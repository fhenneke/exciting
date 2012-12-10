!
!
!
! Copyright (C) 2012 S. Rigamonti and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_writemoke
      Implicit None
Contains
!
!
      Subroutine writemoke (iq, w, moke, fn)
         Use modxs
         Use m_getunit
         Use m_writevars
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq
         Real (8), Intent (In) :: w (:)
         Complex (8), Intent (In) :: moke (:)
         Character (*), Intent (In) :: fn
    ! local variables
         Character (*), Parameter :: thisnam = 'writemoke'
         Integer :: n1 (1), n, iw
         If (any(shape(w) .Ne. shape(moke))) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): input arr&
           &ays have diffenrent shape'
            Call terminate
         End If
         n1 = shape (w)
         n = n1 (1)
         Call getunit (unit1)
         Open (unit1, File=trim(fn), Action='write')
    ! write data to file
         Write (unit1, '(3g18.10)') (w(iw)*escale, moke(iw), iw=1, n)
    ! write relevant parameters to file
         Call writevars (unit1, iq, iq)
         Close (unit1)
      End Subroutine writemoke
!
End Module m_writemoke
