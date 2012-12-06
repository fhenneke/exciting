!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
! !ROUTINE: xslinopt
! !INTERFACE:
Subroutine xslinopt (iq)
! !USES:
      Use modmain
      Use modinput
      Use modxs
      Use modtetra
      Use modmpi
      Use m_genwgrid
      Use m_pade
      Use m_genloss
      Use m_gensigma
      Use m_genmoke
      Use m_gensumrls
      Use m_writeeps
      Use m_writeloss
      Use m_writesigma
      Use m_writemoke
      Use m_writesumrls
      Use m_getunit
      Use m_genfilname
      Implicit None
! !DESCRIPTION:
!   Steps:
!   1- generation of real and imaginary discrete frequency grids
!   2- Reading of inverse dielectric function file for every tensorial component.
!   3- If analytic continuation used, calculate it through pade aproximants.
!   4- symmetrize the macroscopic dielectric function tensor (for optics, i.e. q=0)
!      in order to force the idf to have the exact symmetries of the crystal.
!   5- Generate output file names.
!   6- Calling to subroutines for calculating optical properties.
!   7- Writes calculated properties to files.
!
!   More technical details for each of the above steps:
!   1- Two discrete frequency grids are generated, one of them is on the real axis
!   with a constant small imaginary component added to each point (called wr
!   internally), and the other is purely imaginary (called w internally).
!   2- This subroutine reads the inverse dielectric function from file.
!   The number of tensorial components of the inverse dielectric function (idf)
!   that are read from file,
!   depends on whether q=0 or q!=0. For q=0, all nine components are read
!   (if input%xs%dfoffdiag=True, otherwise only the diagonal components), while
!   for q!=0 only the xx component is read (irrespective of the value of
!   input%xs%dfoffdiag).
!   3- If the read values of the idf come from the evaluation on imaginary
!   frequencies, then these w and wr are used to call the pade soubrutine
!   to generate the values of the idf in the real wr grid by analytic
!   continuation. Otherwise
!   6- Calculated properties are: Loss function, Conductivity, sum rules for optic,
!   moke effect for optic.
! !REVISION HISTORY:
!   Description Dec 2012 (S. Rigamonti)
!EOP
!BOC
  ! arguments
      Integer, Intent (In) :: iq
  ! local variables
      Character (*), Parameter :: thisnam = 'xslinopt'
      Character (256) :: filnam
      Complex (8), Allocatable :: mdf (:), mdf1 (:), mdf2 (:, :, :), w &
     & (:), wr (:), sigma (:), sigma2 (:, :, :),moke(:)
      Real (8), Allocatable :: wplot (:), loss (:)
      Real (8) :: sumrls (3), brd
      Integer :: n, m, recl, iw, nc, oct1, oct2, octl, &
     & octu, optcompt (3)
      Logical :: tq0
      Logical, External :: tqgamma
      tq0 = tqgamma (iq)
  ! number of components (3 for q=0)
      nc = 1
      If (tq0) nc = 3
  ! matrix size for local field effects
      n = ngq (iq)
      Allocate (mdf1(nwdf), mdf2(3, 3, input%xs%energywindow%points), w(nwdf), &
     & wr(input%xs%energywindow%points), wplot(input%xs%energywindow%points), &
     & mdf(input%xs%energywindow%points), loss(input%xs%energywindow%points), &
     & sigma(input%xs%energywindow%points), sigma2(3, 3, input%xs%energywindow%points), &
     & moke(input%xs%energywindow%points))
      mdf2 (:, :, :) = zzero
      sigma2 (:, :, :) = zzero
  ! generate energy grids
      brd = 0.d0
      If (input%xs%tddft%acont) brd = input%xs%broad
  ! w(j) = i*[d*(j-1) + intv(1)], d=[intv(2)-intv(1)]/n, j=1,n, w(n)=intv(2)-d
      Call genwgrid (nwdf, input%xs%energywindow%intv, &
     & input%xs%tddft%acont, 0.d0, w_cmplx=w)
  ! wr(j) = d*(j-1) + intv(1) + i*brd, d=[intv(2)-intv(1)]/n, j=1,n, w(n)=intv(2)-d + i*brd
      Call genwgrid (input%xs%energywindow%points, &
     & input%xs%energywindow%intv, .False., brd, w_cmplx=wr)
  ! w = i*wr + brd
      wplot = dble (wr)
  ! record length
      Inquire (IoLength=Recl) mdf1 (1)
      Call getunit (unit1)
  ! neglect/include local field effects
      Do m = 1, n, Max (n-1, 1)
     ! loop over longitudinal components for optics.
     ! For each tensorial component, read file containing the inverse dielectric function
     ! previously calculated. Fill mdf2 matrix.
         Do oct1 = 1, nc
            If (input%xs%dfoffdiag) Then
               octl = 1
               octu = nc
            Else
               octl = oct1
               octu = oct1
            End If
            Do oct2 = octl, octu
           ! file name for inverse of dielectric function
               Call genfilname (basename='IDF', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctype=input%xs%tddft%fxctypenumber, tq0=tq0, &
              & oc1=oct1, oc2=oct2, iqmt=iq, filnam=filnam)
           ! read macroscopic dielectric function (original frequencies)
               Open (unit1, File=trim(filnam), Form='unformatted', &
              & Action='read', Status='old', Access='direct', &
              & Recl=Recl)
               Do iw = 1, nwdf
                  Read (unit1, Rec=iw) mdf1 (iw)
               End Do
               Close (unit1)
           ! analytic continuation
               If (input%xs%tddft%acont) Then
                  Call pade (input%xs%energywindow%points, wr, nwdf, w, &
                 & mdf1, mdf)
               Else
                  mdf (:) = mdf1 (:)
               End If
               mdf2 (oct1, oct2, :) = mdf (:)
            End Do
         End Do
     ! loop again over longitudinal components for optics.
         Do oct1 = 1, nc
            If (input%xs%dfoffdiag) Then
               octl = 1
               octu = nc
            Else
               octl = oct1
               octu = oct1
            End If
            Do oct2 = octl, octu
               optcompt (:) = (/ oct1, oct2, 0 /)
           ! symmetrize the macroscopic dielectric function tensor.
           ! This forces the dielectric tensor to have the symmetries of the crystal.
           ! Note: This symmetrization is not related to the one described in
           ! chapter 8.3 of S.Sagmeister thesis.
               if (tq0) Call symt2app (oct1, oct2, nwdf, symt2, mdf2, mdf)
           ! file names for spectra
               Call genfilname (basename='EPSILON', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctypestr=input%xs%tddft%fxctype, tq0=tq0, &
              & oc1=oct1, oc2=oct2, iqmt=iq, filnam=fneps)
               Call genfilname (basename='LOSS', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctypestr=input%xs%tddft%fxctype, tq0=tq0, &
              & oc1=oct1, oc2=oct2, iqmt=iq, filnam=fnloss)
               Call genfilname (basename='SIGMA', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctypestr=input%xs%tddft%fxctype, tq0=tq0, &
              & oc1=oct1, oc2=oct2, iqmt=iq, filnam=fnsigma)
               if (tq0) Call genfilname (basename='SUMRULES', asc=.False., &
              & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
              & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
              & fxctypestr=input%xs%tddft%fxctype, tq0=tq0, &
              & oc1=oct1, oc2=oct2, iqmt=iq, filnam=fnsumrules)

           ! generate optical functions
               Call genloss (mdf, loss)
               Call gensigma (wplot, mdf, optcompt, sigma)
               sigma2 (oct1, oct2, :) = sigma (:)
               Call gensumrls (wplot, mdf, sumrls)
           ! write optical functions to file
               Call writeeps (iq, oct1, oct2, wplot, mdf, trim(fneps))
               Call writeloss (iq, wplot, loss, trim(fnloss))
               Call writesigma (iq, wplot, sigma, trim(fnsigma))
               if (tq0) Call writesumrls (iq, sumrls, trim(fnsumrules))
           ! end loop over optical components
            End Do
         End Do
        If (tq0) Then
            Call genfilname (basename='MOKE', asc=.False., &
            & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
            & input%xs%tddft%aresdf, tord=input%xs%tddft%torddf, nlf=(m == 1), &
            & fxctypestr=input%xs%tddft%fxctype, tq0=tq0, iqmt=iq, filnam=fnmoke)
            Call genmoke (input%xs%energywindow%points,wplot,sigma2, moke)
            Call writemoke (iq, wplot, moke, trim(fnmoke))
        End If
      End Do ! m
  ! deallocate
      Deallocate (mdf, mdf1, mdf2, w, wr, wplot, loss, sigma, sigma2, moke)
End Subroutine xslinopt
