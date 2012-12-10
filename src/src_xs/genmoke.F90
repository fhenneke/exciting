!
!
!
! Copyright (C) 2012 S. Rigamonti and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_genmoke
      Implicit None
Contains
!
!
!BOP
! !ROUTINE: genmoke
! !INTERFACE:
        Subroutine genmoke(n, w, sigma, moke)
! !USES:
         Use modxs
         Use mod_constants, Only: pi, zi
! !DESCRIPTION:
!   Calculates the complex kerr parameters (or angle) for polar geometry.
!   The calculation is done according to Eq. 1 of Ref. PRB 45, 10924 (1992).
! !REVISION HISTORY:
!   Created Dec 2012 (S. Rigamonti)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: n
         Real (8), Intent (In) :: w (:)
         Complex (8), Intent (In) :: sigma (:,:,:)
         Complex (8), Intent (Out) :: moke (:)
    ! local variables
         Character (*), Parameter :: thisnam = 'genmoke'
         Real (8), Allocatable :: ones (:)

         Allocate (ones(n))
         ones(:)=1.d0

         moke (:) = sqrt( ones(:) + zi * 4.d0*pi * sigma(1,1,:) / w (:) )

         moke (:) = - sigma(1,2,:) / ( sigma(1,1,:)*moke (:) )

         Return
    End Subroutine genmoke

End Module m_genmoke
