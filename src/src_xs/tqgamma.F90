!
!
!
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: tqgamma
! !INTERFACE:
Logical Function tqgamma (iq)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   iq          : index of qpoint (in,integer)
!   tqgamma     : S (solution) matrix (out,boolean)!
!DESCRIPTION:
! Checks whether the q point corresponding to the index iq is the
! Gamma point or not. If this is the case, the tqgamma variable is
! set to true.
! !REVISION HISTORY:
!   Description Dec 2012 (S. Rigamonti)
!EOP
!BOC
  ! arguments
      Integer, Intent (In) :: iq
  ! local variables
      Real (8) :: epsg = 1.d-12
      tqgamma = .False.
      If (sum(Abs(vqc(:, iq))) .Lt. epsg) tqgamma = .True.
End Function tqgamma
