!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: df
! !INTERFACE:
!
!
Module m_gensigma
      Implicit None
Contains
!
!
      Subroutine gensigma (w, eps, oc, sigma)
         Use mod_constants, Only: pi, zi
         Use modxs
         Use modinput
! !DESCRIPTION:
!   Calculation of conductivity tensor (sigma) follows Eq.1 from
!   PRB 86, 125139 (2012)
! !REVISION HISTORY:
!   Added description Nov 2012 (S. Rigamonti)
!EOP
!BOC
         Implicit None
    ! arguments
         Real (8), Intent (In) :: w (:)
         Complex (8), Intent (In) :: eps (:)
         Integer, Intent (In) :: oc (3)
         Complex (8), Intent (Out) :: sigma (:)
    ! local variables
         Character (*), Parameter :: thisnam = 'gensigma'
         Real (8) :: delt
         Real (8) :: wp,gd
         wp = 0.275619816d0
         gd = 0.02d0
         If (any(shape(eps) .Ne. shape(sigma))) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): input and&
           & output arrays have diffenrent shape'
            Call terminate
         End If
    ! optical conductivity
         delt = 0.0d0
         If (oc(1) .Eq. oc(2)) delt = 1.d0
         sigma (:) = -zi * w(:) * (eps(:) - delt) / (4.d0*pi)
!         sigma (:) = aimag (eps(:)) * w (:) / (4.d0*pi)
!         sigma (:) = sigma (:) + zi * &
!        & (-(dble(eps(:))-delt)*w(:)/(4.d0*pi))
         If (input%xs%tddft%intraband .And. (oc(1) .Eq. oc(2))) Then
            write(*,*) "calculating drude term in diagonal component of sigma"
            sigma (:) = sigma (:) + wp*wp/(4.d0*pi*(gd - zi*w(:)))
         End If
      End Subroutine gensigma
!
End Module m_gensigma
