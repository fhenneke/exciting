!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine putoccsv (ik, occsvp)
      Use modmain
      Use modmpi
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Real (8), Intent (In) :: occsvp (nstsv)
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
!
! local variables
      Integer :: recl, koffset
! find the record length
      Inquire (IoLength=Recl) vkl (:, ik), nstsv, occsvp
!$OMP CRITICAL
      filetag = 'OCCSV'
      If (splittfile .Or. (rank .Eq. 0)) Then
         Open (70, File=outfilenamestring(filetag, ik), Action='WRITE', &
        & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
         If (splittfile) Then
            koffset = ik - firstk (procofk(ik)) + 1
         Else
            koffset = ik
         End If
         Write (70, Rec=koffset) vkl (:, ik), nstsv, occsvp
         Close (70)
!
      End If
!$OMP END CRITICAL
      Return
End Subroutine
