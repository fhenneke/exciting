#!/usr/bin/env python
# -*- coding: utf-8 -*-
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%% -------------------------------------------- OPTIMIZE-analyze.py -------------------------------------------- %%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#
# AUTHOR:
#    Rostam Golesorkhtabar  
#    r.golesorkhtabar@gmail.com
# 
# DATE:
#    Wed Aug 01 00:00:00 2012
#
# SYNTAX:
#    python OPTIMIZE-analyze.py
#           OPTIMIZE-analyze.py
# 
# EXPLANATION:
#
#______________________________________________________________________________________________________________________

from   pylab import *
import os
import sys
import glob
import copy
import math
import os.path
import numpy as np
from   lxml  import etree as ET
import matplotlib.pyplot as plt
import pylab             as pyl
from   scipy.optimize import fmin_powell

#%%%%%%%%--- CONSTANTS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
_e   =  1.602176565e-19            # elementary charge
Bohr =  5.291772086e-11            # a.u. to meter
Ry2eV= 13.605698066                # Ry to eV
cnvrtr = (_e*Ry2eV)/(1e9*Bohr**3)  # Ry/[a.u.]^3 to GPa
#----------------------------------------------------------------------------------------------------------------------

#%%%--- SUBROUTINS AND FUNCTIONS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def M(p0, V):
    """ Murnaghan """
    E0, V0, B0, Bp = p0
    E  = E0 + (B0*V/Bp*(1/(Bp-1)*(V0/V)**Bp +1)-B0*V0/(Bp-1))
    return E
#----------------------------------------------------------------------------------------------------------------------

def BM(p0, V):
    """ Birch-Murnaghan """
    E0, V0, B0, Bp = p0
    E = E0+(9.*B0*V0/16)*(((((V0/V)**(2./3))-1.)**3.)*Bp + ((((V0/V)**(2/3.))-1.)**2.)*(6.-4.*((V0/V)**(2./3.))))
    return E
#----------------------------------------------------------------------------------------------------------------------

def snr(p0, v, e):
    """ Squared norm of residue vector calculation """
    if (fit == 'M'): 
        E_M = M(p0, v)
        return np.sum((e - E_M)**2.)
    else:
        E_BM = BM(p0, v)
        return np.sum((e - E_BM)**2.)
#----------------------------------------------------------------------------------------------------------------------

def sortlist(lst1, lst2):
    temp = copy.copy(lst1)

    lst3 = []
    lst4 = []

    temp.sort()

    for i in range(len(lst1)):
        lst3.append(lst1[lst1.index(temp[i])])
        lst4.append(lst2[lst1.index(temp[i])])

    return lst3, lst4
#______________________________________________________________________________________________________________________

#%%%--- Reading the INFO file ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INFO_file = str(glob.glob('INFO_*')[0])

INFO=open(INFO_file,'r')
mod =INFO_file[5:]

l1  = INFO.readline()
l2  = INFO.readline()

l3  = INFO.readline()
mdr = float(l3.split()[-1])

l4  = INFO.readline()
NoP = int(l4.split()[-1])

INFO.close()
#----------------------------------------------------------------------------------------------------------------------

#%%%--- plot definitions    ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

emt='%11.7f'
fmt='%8.4f'
bmt='%7.3f'
pmt='%6.1f'

params = {'axes.linewidth'        : 2.,
          'figure.subplot.bottom' : 0.14,
          'figure.subplot.right'  : 0.93, 
          'figure.subplot.left'   : 0.20,
          'xtick.major.pad'       : 8, 
          }

plt.rcParams.update(params)

fig = plt.figure()
ax = fig.add_subplot(111)

#%%%--- Reading the energies ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (mod == 'VOL'):

    volume = []
    energy = []

    vollist= glob.glob('vol_??')
    for vol_num in vollist:
        os.chdir(vol_num)

        if (os.path.exists('INFO.OUT') == False):
            print '\n     ... Oops NOTICE: There is NO "INFO.OUT" file !?!?!?    \n'

        for line in open('INFO.OUT','r'):
            if (line.find('Unit cell volume')>=0): 
                vol = float(line.split()[-1])
                break
        volume.append(vol)

        if (os.path.exists('TOTENERGY.OUT') == False):
            print '\n     ... Oops NOTICE: There is NO "TOTENERGY.OUT" file !?!?!?    \n'

        e_file = open('TOTENERGY.OUT', 'r')
        enrgis = e_file.readlines()
        ene = float(enrgis[-1])
        e_file.close()
        energy.append(ene)

        os.chdir('../')

    volume, energy = sortlist(volume, energy)

    fvol = open('energy-vs-volume.dat', 'w')
    for i in range(len(energy)):
        print >>fvol, volume[i],'   ', energy[i]
    fvol.close()

    data = np.loadtxt('energy-vs-volume.dat')
    vi, ei = data.T
    if (len(ei) < 3): sys.exit('\n     ... Oops ERROR: EOS fit needs at least 3 points.    \n')
    ei = ei*2.

    #%%%%%%%%--- Reading the fit type ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fit = raw_input('\n>>>> Murnaghan or Birch-Murnaghan fit: [M/B] ').upper()
    if (fit != 'B' and fit != 'M'): sys.exit("\n.... Oops ERROR: Choose 'B' or 'M' \n")

    #%%%%%%%%--- FIT CALCULATIONS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a2, a1, a0 = np.polyfit(vi, ei, 2)
    V0 = -a1/(2.*a2)
    E0 = a2*V0**2. + a1*V0 + a0
    B0 = a2*V0
    Bp = 2.
     
    p0 = [E0, V0, B0, Bp]

    viei = sorted([zip(vi, ei)])
    v, e = np.array(viei).T

    p1, fopt, direc, n_iter, n_funcalls, warnflag = fmin_powell(snr, p0, args=(v, e), full_output=True, disp=0)
    E0, V0, B0, Bp = p1

    E0_dim = E0/2.
    B0_GPa = B0*cnvrtr
    fopt = fopt/(4.)
    fopt = sqrt(fopt)
        
    print\
    '\n Log [Final residue in Hartree]:',pmt%(log10(fopt)),'\n'\
    '\n === Final parameters ==='      \
    '\n E0 = ',emt%(E0_dim), '  [Ha]'\
    '\n V0 = ',fmt%(V0),' [a.u.^3]'\
    '\n B0 = ',bmt%(B0_GPa),'     [GPa]'\
    "\n B' = ",bmt%(Bp),\
    "\n ========================\n"
    
    vmn = min(np.min(vi), V0)
    vmx = max(np.max(vi), V0)
    dv  = vmx - vmn
    vfit = np.linspace(vmn-(0.1*dv), vmx+(0.1*dv), 1000)

    if (fit=='M'): efit = M(p1, vfit)
    if (fit=='B'): efit =BM(p1, vfit)

    #%%%--- Writing the 'vol-optimized.xml' file ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    INOBJ= open(mod.lower()+'-xml/input.xml', 'r')
    doc  = ET.parse(INOBJ)
    root = doc.getroot()

    scale = map(float,doc.xpath('/input/structure/crystal/@scale'))
    if (scale==[]):
        sys.exit('\n     ... Oops ERROR: There is NO scale in '+ INF +' file !?!?!?    \n')

    stretchstr = doc.xpath('/input/structure/crystal/@stretch')
    if (stretchstr==[]):
        stretch=[1.,1.,1.]
    else:
        stretch=np.array(map(float,stretchstr[0].split()))

    basevectsn = doc.xpath('//basevect/text()')
    bv = []
    for basevect in basevectsn:
        bv.append(map(float,basevect.split()))

    M_old= np.array(bv)
    D    = np.linalg.det(M_old)
    V0_in= abs(stretch[0]*stretch[1]*stretch[2]*scale[0]**3*D)
 
    crystal = doc.xpath('//crystal')
    new_scale = scale[0]*(V0/V0_in)**(1./3.)
    crystal[0].set("scale",str(round(new_scale,10)))

    # Writing the structure file --------------------------------------------------------------------------------------
    OUTOBJ   = open('vol-optimized.xml', 'w')
    OUTOBJ.write(ET.tostring(root, method         ='xml',
                                   pretty_print   =True ,
                                   xml_declaration=True ,
                                   encoding       ='UTF-8'))
    OUTOBJ.close()

    #%%%%%%%%--- PLOT PREPARATION ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    if (fit=='M'): fit_label = 'Murnaghan fit'
    if (fit=='B'): fit_label = 'Birch-Murnaghan fit'
   
    xlabel = u'Volume [Bohr\u00B3]'
    ylabel = 'Energy [Ha]'

    plt.text(0.35,0.75, 'E$_{min}$ = '+str(round(E0_dim,7))+' [Ha]',      transform = ax.transAxes)
    plt.text(0.35,0.70, 'V$_{min}$ = '+str(round(V0,4))+u' [Bohr\u00B3]', transform = ax.transAxes)
    plt.text(0.35,0.65, 'B$_0$ = '+str(round(B0_GPa,3))+' [GPa]',     transform = ax.transAxes)
    plt.text(0.35,0.60, 'B$^\prime$  = '+str(round(Bp,3))       ,     transform = ax.transAxes)

    xx = [] ; xx = vfit
    yy = [] ; yy = efit
    x0 = [] ; x0 = vi
    y0 = [] ; y0 = ei

#----------------------------------------------------------------------------------------------------------------------
if (mod != 'VOL'):
    print
    fee  = open('energy-vs-strain.dat', 'w')

    for i in range(1, NoP+1):
        if (i<10):
            dir_num = mod.lower() +'_0'+str(i)
        else:
            dir_num = mod.lower() +'_' +str(i)

        if (os.path.exists(dir_num) == False):
            print '\n     ... Oops NOTICE: There is NO '+ dir_num +' directory !?!?!?    \n'
            break
    
        os.chdir(dir_num)

        if (os.path.exists('TOTENERGY.OUT') == False):
            print '\n     ... Oops NOTICE: There is NO "TOTENERGY.OUT" file !?!?!?    \n'

        s = i-(NoP+1)/2   
        r = 2*mdr*s/(NoP-1)
        if (s==0): r=0.00001

        if (r>0):
            strain ='+'+str(round(r,10))
        else:
            strain = str(round(r,10))

        if (os.path.exists('TOTENERGY.OUT')):
            e_file = open('TOTENERGY.OUT', 'r')
            enrgis = e_file.readlines()
            energy = float(enrgis[-1])
            e_file.close()

        print >>fee, strain,'   ', energy
        os.chdir('../')

    fee.close()

    data = np.loadtxt('energy-vs-strain.dat')
    si, ei = data.T
    vs = sorted([zip(si, ei)])
    s, e = np.array(vs).T
    if (len(e) < 5): sys.exit('\n     ... Oops ERROR: 4th order polynomial fit needs at least 5 points.\n')

    #%%%%%%%%--- FIT CALCULATIONS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coeffitions = np.polyfit(si, ei, 4)
    f4 = np.poly1d(coeffitions)

    sfit = np.linspace(mdr*-1.2, mdr*1.2, 1000)
    efit = f4(sfit)
    s_min= sfit[efit.argmin()]
    e_min= efit[efit.argmin()]

    #%%%--- Writing the 'mod-optimized.xml' file ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    INOBJ= open(mod.lower()+'-xml/input.xml', 'r')
    doc  = ET.parse(INOBJ)
    root = doc.getroot()

    stretchstr = doc.xpath('/input/structure/crystal/@stretch')
    if (stretchstr==[]):
        stretch=[1.,1.,1.]
    else:
        stretch=np.array(map(float,stretchstr[0].split()))

    basevectsn = doc.xpath('//basevect/text()')
    bv = []
    for basevect in basevectsn:
        bv.append(map(float,basevect.split()))

    M_old= np.array(bv)

    def_matrix={\
    'BOA'  :[[(1+s_min)**-.5, 0.          , 0.           ],
             [ 0.           , 1.+s_min    , 0.           ],
             [ 0.           , 0.          ,(1+s_min)**-.5]],\

    'COA'  :[[(1+s_min)**-.5, 0.            , 0.         ],
             [ 0.           , (1+s_min)**-.5, 0.         ],
             [ 0.           , 0.            , 1.+s_min   ]],\

    'ALPHA':[[1./(1-s_min**2), 0.           , 0.         ],
             [ 0.            , 1.           ,s_min       ],
             [ 0.            ,s_min         , 1.         ]],\

    'BETA' :[[ 1.          , 0.             ,s_min       ],
             [ 0.          , 1./(1-s_min**2), 0.         ],
             [s_min        , 0.             , 1.         ]],\

    'GAMMA':[[ 1.          ,s_min       , 0.             ],
             [s_min        , 1.         , 0.             ],
             [ 0.          , 0.         , 1./(1-s_min**2)]]}

    M_smn = np.array(def_matrix[mod])
    M_new = np.dot(M_old, M_smn)

    bsvct = doc.xpath('//crystal/basevect')
    for j in range(3):
        bdummy = '%22.16f'%(M_new[j,0]) + '%22.16f'%(M_new[j,1]) + '%22.16f'%(M_new[j,2])+' ' 
        bsvct[j].text = bdummy

    #---Writing the structure file-------------------------------------------------------------------------------------
    OUTOBJ   = open(mod.lower()+'-optimized.xml', 'w')
    OUTOBJ.write(ET.tostring(root, method         ='xml',
                                   pretty_print   =True ,
                                   xml_declaration=True ,
                                   encoding       ='UTF-8'))
    OUTOBJ.close()

    #%%%%%%%%--- PLOT PREPARATION ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xlabel      = 'Physical strain $\epsilon$'
    ylabel      = 'Energy [Ha]'
    fit_label   = '4th order polynomial fit'

    xx = [] ; xx = sfit
    yy = [] ; yy = efit
    x0 = [] ; x0 = si
    y0 = [] ; y0 = ei
    
    plt.text(0.35,0.75, 'E$_{min}$ = '+str(round(e_min, 7))+' [Ha]', transform = ax.transAxes)
    plt.text(0.35,0.70, '$\epsilon_{min}$  = '+str(round(s_min,5)), transform = ax.transAxes)
    
#----------------------------------------------------------------------------------------------------------------------

print " Optimized structure saved to file:",mod.lower()+"-optimized.xml\n"

ax.set_xlabel(xlabel, fontsize = 18)
ax.set_ylabel(ylabel, fontsize = 18)
ax.plot(xx, yy, 'k', color='red'  , linewidth =2, label=fit_label)
ax.plot(x0, y0, 'o', color='green', markersize=8, markeredgecolor='black',markeredgewidth = 1,label='Calculations')
ax.legend(numpoints=1,loc=9)

for label in ax.xaxis.get_ticklabels(): label.set_fontsize(15)
for label in ax.yaxis.get_ticklabels(): label.set_fontsize(15)
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)

pyl.grid(True)   

ax.xaxis.set_major_locator(MaxNLocator(7))

dyy = (max(yy)-min(yy))/15
ax.set_ylim(min(yy)-dyy,max(yy)+dyy)
dxx = (max(xx)-min(xx))/18
ax.set_xlim(min(xx)-dxx,max(xx)+dxx)

plt.savefig(mod.lower()+'.png', orientation='portrait',format='png',dpi=80)
plt.savefig(mod.lower()+'.ps', orientation='portrait',format='eps')
plt.show()
#----------------------------------------------------------------------------------------------------------------------
