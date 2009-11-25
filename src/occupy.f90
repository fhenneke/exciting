!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: occupy
! !INTERFACE:
!
!
Subroutine occupy
! !USES:
      Use modinput
      Use modmain
#ifdef TETRAOCC_DOESNTWORK
!<rga>
      Use modtetra
!</rga>
#endif
! !DESCRIPTION:
!   Finds the Fermi energy and sets the occupation numbers for the
!   second-variational states using the routine {\tt fermi}.
!
! !REVISION HISTORY:
!   Created February 2004 (JKD)
!   Modifiactions for tetrahedron method, November 2007 (RGA alias
!     Ricardo Gomez-Abal)
!EOP
!BOC
      Implicit None
! local variables
      Integer, Parameter :: maxit = 1000
      Integer :: ik, ist, it
      Real (8) :: e0, e1, chg, x, t1
! external functions
      Real (8) :: sdelta, stheta
      External sdelta, stheta
! find minimum and maximum eigenvalues
      e0 = evalsv (1, 1)
      e1 = e0
      Do ik = 1, nkpt
         Do ist = 1, nstsv
            e0 = Min (e0, evalsv(ist, ik))
            e1 = Max (e1, evalsv(ist, ik))
         End Do
      End Do
      If (e0 .Lt. input%groundstate%evalmin) Then
         Write (*,*)
         Write (*, '("Warning(occupy): valence eigenvalues below evalmi&
        &n for s.c. loop ", I5)') iscl !"
      End If
#ifdef TETRAOCC_DOESNTWORK
!<rga>
      If ( .Not. istetraocc()) Then
!</rga>
#endif
         t1 = 1.d0 / input%groundstate%swidth
! determine the Fermi energy using the bisection method
         Do it = 1, maxit
            efermi = 0.5d0 * (e0+e1)
            chg = 0.d0
            Do ik = 1, nkpt
               Do ist = 1, nstsv
                  If (evalsv(ist, ik) .Gt. input%groundstate%evalmin) &
                 & Then
                     x = (efermi-evalsv(ist, ik)) * t1
                     occsv (ist, ik) = occmax * stheta &
                    & (input%groundstate%stypenumber, x)
                     chg = chg + wkpt (ik) * occsv (ist, ik)
                  Else
                     occsv (ist, ik) = 0.d0
                  End If
               End Do
            End Do
            If (chg .Lt. chgval) Then
               e0 = efermi
            Else
               e1 = efermi
            End If
            If ((e1-e0) .Lt. input%groundstate%epsocc) Go To 10
         End Do
         Write (*,*)
         Write (*, '("Error(occupy): could not find Fermi energy")')
         Write (*,*)
         Stop
10       Continue
! find the density of states at the Fermi surface in units of
! states/Hartree/spin/unit cell
         fermidos = 0.d0
         Do ik = 1, nkpt
            Do ist = 1, nstsv
               If (evalsv(ist, ik) .Gt. input%groundstate%evalmin) Then
                  x = (evalsv(ist, ik)-efermi) * t1
                  fermidos = fermidos + wkpt (ik) * sdelta &
                 & (input%groundstate%stypenumber, x) * t1
               End If
            End Do
            If (occsv(nstsv, ik) .Gt. input%groundstate%epsocc) Then
               Write (*,*)
               Write (*, '("Warning(occupy): not enough empty states fo&
              &r k-point ", I6)') ik
            End If
         End Do
         fermidos = fermidos * occmax
#ifdef TETRAOCC_DOESNTWORK
!<rga>
      Else
  ! calculate the Fermi energy and the density of states at the Fermi energy
         Call fermitetifc (nkpt, nstsv, evalsv, chgval, &
        & associated(input%groundstate%spin), efermi, fermidos)
         Call tetiwifc (nkpt, nstsv, evalsv, efermi, occsv)
         Do ik = 1, nkpt
    ! The "occsv" variable returned from "tetiw" already contains the
    ! weight "wkpt" and does not account for spin degeneracy - rescaling is
    ! necessary (S. Sagmeister).
            Do ist = 1, nstsv
               occsv (ist, ik) = (occmax/wkpt(ik)) * occsv (ist, ik)
            End Do
         End Do
      End If
!</rga>
#endif
      Return
End Subroutine
!EOC
