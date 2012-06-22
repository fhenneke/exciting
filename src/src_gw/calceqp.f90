!BOP
!
! !ROUTINE: calceqp
!
! !INTERFACE:
      subroutine calceqp

! !DESCRIPTION:
! 
! Given the matrix elements $\langle
! \Psi_{n\vec{k}}|\Sigma(\vec{k},\omega)|\Psi_{n\vec{k}}\rangle$, $\langle
! \Psi_{n\vec{k}}|V^{xc}|\Psi_{n\vec{k}}\rangle$ and
! $\varepsilon^{DFT}_{n\vec{k}}$. this subroutine calculates the
! quasi-particle energies $\varepsilon^{qp}_{n\vec{k}}$
!
! !USES:
!
      use modmain
      use modgw    
       
! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: ie   !(Counter) Runs over bands
      integer(4) :: ikp  !(Counter) Runs over k-points
      integer(4) :: npar ! Number of parameters of the analitic cont. of
!                          the selfenergy
      integer(4) :: it
      integer(4) :: ierr
      
      real(8) :: enk,snk,vxcnk,znk
      real(8) :: delta 

      real(8) :: eferqp,es
      real(8) :: egap,egap0,egapold
      
      complex(8) :: ein 
      complex(8) :: dsig
      complex(8) :: sigma
      
      complex(8), allocatable :: sigc(:,:)
      real(8), allocatable    :: znorm(:,:)

      complex(8), allocatable :: a(:),sc(:),poles(:) ! Parameters for AC

! !DEFINED PARAMETERS:

      real(8), parameter :: etol=1.0d-6
      integer, parameter :: nitmax=100

! A critical parameters used in this subroutine is iopes, 
! which control how Fermi energy shift is treated 
!     iopes =    
!          0 -- perturbative G0W0 without energy shift
!          1 -- perturbative G0W0 with energy shift
!          2 -- iterative G0W0 with energy shift
!          3 -- iterative G0W0 without energy shift
 
! !EXTERNAL ROUTINES: 

      external fermi
      external writeqp


! !REVISION HISTORY:
!
! Created: 16.08.05 by RGA
! Revisited June 2011 by DIN
!
!EOP
!BOC

      call boxmsg(fgw,'-','QP-energy output')
!
!     Allocate the array for the quasi-particle energies 
! 
      allocate(eqp(ibgw:nbgw,nkpt))
      eqp(ibgw:nbgw,1:nkpt)=evaldft(ibgw:nbgw,1:nkpt)
      
      allocate(sigc(ibgw:nbgw,nkpt))
      allocate(znorm(ibgw:nbgw,nkpt))
!
!     Calculate the number of parameters of the analitic function for the selfenergy
!      
      npar=2*npol
      allocate(a(npar),sc(nomeg),poles(npar))
      allocate(sacpar(npar,ibgw:nbgw,nkpt))

!----------------------------------
!     Start iterative procedure
!----------------------------------      
      eferqp=0.0
      egap=0.0
      ierr=0

      do it = 0, nitmax
      
         do ikp = 1, nkpt
           do ie = ibgw, nbgw
           
              enk=evaldft(ie,ikp)
!          
!            Set up the analytic continuation parameters 
!
              if(it.eq.0)then
                 sc=selfec(ie,ikp,1:nomeg)
                 call setsac(iopac,nomeg,npar,enk,sc,freqs,a,poles)
                 sacpar(:,ie,ikp)=a
              else
                 a=sacpar(:,ie,ikp)
              end if
              
              if (iopes.eq.2) then 
                ein=cmplx(eqp(ie,ikp)-es,0.0d0,8)
              else if (iopes.eq.3) then 
                ein=cmplx(eqp(ie,ikp),0.0d0,8)
              else 
                ein=cmplx(enk,0.0d0,8)
              endif 
!
!             Perform AC
!
              ein=cmplx(enk,0.0d0,8)
              call getsac(iopac,nomeg,npar,enk,ein,freqs,a,sigma,dsig)

!             Set the correlation selfenergy
              sigc(ie,ikp)=sigma

!             Set the normalization factor              
              znk=1.0d0/(1.0d0-real(dsig))
              if(znk.gt.1.d0 .or. znk.lt.0.d0) then
               write(fgw,*)'WARNING(calceqp): nonphysical Znk for ikp=', ikp
               znk=0.8 
              endif
              znorm(ie,ikp)=znk
!            
!             Calculate the new quasi-particle energy 
!
              snk=real(selfex(ie,ikp))+real(sigc(ie,ikp))
              vxcnk=real(vxcnn(ie,ikp))
              
              select case(iopes) 
              case(0)
                delta=znk*(snk-vxcnk)
              case(1) 
                delta=znk*(snk-vxcnk)+(1.0d0-znk)*es
              case (2) 
                if(it.eq.0) then
                  delta=znk*(snk-vxcnk)
                else
                  delta=snk-vxcnk
                endif
              case(3) 
                delta=snk-vxcnk
              end select 

              eqp(ie,ikp)=enk+delta
          
           enddo ! ie
         enddo ! ikp

         call fermi(nkpt,nbandsgw,eqp,nirtet,tndi,wirtet,tvol, &
        &    nvelgw,.false.,eferqp,egap)

         if(it.eq.0) then 
            egap0=egap
            write(fgw,8) 
         endif 
       
         if(egap.gt.0.0d0) then                                                
            write(fgw,9) it,eferqp,es,egap*hev,(eqp(minunoband,1)-eqp(maxoccband,1))*hev
         else                                                                 
            write(fgw,*) 'WARNING(calceqp):!!! metallic, DOS at Fermi level', -egap                
         endif 
    
         if(it.ne.0)then 
           if(abs(egap-egap0).gt.0.2d0) then
             write(fgw,*) 'WARNING(calceqp): --- Band gap deviates from initial GW gap more than 0.2 Ha'
             ierr=1
             exit 
           end if
         end if

         if( (iopes.eq.0).or. &
        &   ((abs(es-eferqp).lt.etol).and.(abs(egap-egapold).le.etol)) ) exit

!        Perform next iteration
         es=eferqp
         egapold=egap
         
      enddo ! it
      
      if(it.gt.nitmax)ierr=1
      if(ierr.ne.0) then 
         write(fgw,*) 'WARNING(calceqp): --- Failed to converge!!!'
      endif 
!      
!     Write quasi-particle energies to disk
!      
      call writeqp(sigc,znorm)
      
      deallocate(a)
      deallocate(znorm)
      deallocate(sigc)

!     KS band structure
      call bandanaly(ibgw,nbgw,nkpt,vkl,evaldft(ibgw:nbgw,:),efermi,"KS",fgw)

!     QP band structure
      call bandanaly(ibgw,nbgw,nkpt,vkl,eqp(ibgw:nbgw,:),eferqp,"G0W0",fgw)

   8  format( ' #iter',7x,"Ef_QP",10x,"es",9x,"Eg/eV",2x,'Eg(k=0)/eV')
   9  format( i6,4f12.6)

      return
      end subroutine calceqp

!EOC        
          