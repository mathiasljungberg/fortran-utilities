!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 
!    02110-1301  USA
!

LOGICAL FUNCTION match_wild(pattern, string)
   ! compare given string for match to pattern which may
   ! contain wildcard characters:
   ! "?" matching any one character, and
   ! "*" matching any zero or more characters.
   ! Both strings may have trailing spaces which are ignored.
   ! Authors: Clive Page, userid: cgp  domain: le.ac.uk, 2003 (original code)
   !          Rolf Sander, 2005 (bug fixes and pattern preprocessing)
   ! Minor bug fixed by Clive Page, 2005 Nov 29, bad comment fixed 2005 Dec 2.

!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 
!    02110-1301  USA
!
   IMPLICIT NONE

   CHARACTER(LEN=*), INTENT(IN) :: pattern ! pattern may contain * and ?
   CHARACTER(LEN=*), INTENT(IN) :: string  ! string to be compared
   INTEGER :: lenp, lenp2, lens, n, p2, p, s
   INTEGER :: n_question, n_asterisk

   CHARACTER(LEN=LEN(pattern)) :: pattern2

   lens = LEN_TRIM(string)
   lenp = LEN_TRIM(pattern)

! If the pattern is empty, always return true
   IF (lenp == 0) THEN
     match_wild = .TRUE.
     RETURN
   ENDIF

! The pattern must be preprocessed. All consecutive occurences of
! one or more question marks ('?') and asterisks ('*') are sorted and
! compressed. The result is stored in pattern2.

   pattern2(:)=''
   p  = 1 ! current position in pattern
   p2 = 1 ! current position in pattern2
   DO
     IF ((pattern(p:p) == '?').OR.(pattern(p:p) == '*')) THEN
! a special character was found in the pattern
       n_question = 0
       n_asterisk = 0
       DO WHILE (p <= lenp)
         ! count the consecutive question marks and asterisks
         IF ((pattern(p:p) /= '?').AND.(pattern(p:p) /= '*')) EXIT
         IF (pattern(p:p) == '?') n_question = n_question + 1
         IF (pattern(p:p) == '*') n_asterisk = n_asterisk + 1
         p = p + 1
       ENDDO
       IF (n_question>0) THEN ! first, all the question marks
         pattern2(p2:p2+n_question-1) = REPEAT('?',n_question)
         p2 = p2 + n_question
       ENDIF
       IF (n_asterisk>0) THEN ! next, the asterisk (only one!)
         pattern2(p2:p2) = '*'
         p2 = p2 + 1
       ENDIF
     ELSE
! just a normal character
       pattern2(p2:p2) = pattern(p:p)
       p2 = p2 + 1
       p = p + 1
     ENDIF
     IF (p > lenp) EXIT
   ENDDO
!!   lenp2 = p2 - 1
   lenp2 = len_trim(pattern2)

! The modified wildcard in pattern2 is compared to the string:

   p2 = 1
   s = 1
   match_wild = .FALSE.
   DO
     IF (pattern2(p2:p2) == '?') THEN
! accept any char in string
       p2 = p2 + 1
       s = s + 1
     ELSEIF (pattern2(p2:p2) == "*") THEN
       p2 = p2 + 1
       IF (p2 > lenp2) THEN
! anything goes in rest of string
         match_wild = .TRUE.
         EXIT ! .TRUE.
       ELSE
! search string for char at p2
         n = INDEX(string(s:), pattern2(p2:p2))
         IF (n == 0) EXIT  ! .FALSE.
         s = n + s - 1
       ENDIF
     ELSEIF (pattern2(p2:p2) == string(s:s)) THEN
! single char match
       p2 = p2 + 1
       s = s + 1
     ELSE
       ! non-match
       EXIT ! .FALSE.
     ENDIF
     IF (p2 > lenp2 .AND. s > lens) THEN
! end of both pattern2 and string
       match_wild = .TRUE.
       EXIT ! .TRUE.
     ENDIF

!!     IF (s > lens .AND. (pattern2(p2:p2) == "*") .AND. p2 == lenp2) THEN
!! above line buggy since p2 can be beyond end of string pattern2 by this point. CGP

       IF (s > lens .AND. p2 == lenp) THEN
         IF(pattern2(p2:p2) == "*") THEN
! "*" at end of pattern2 represents an empty string
           match_wild = .TRUE.
           EXIT
         END IF
     ENDIF
     IF (p2 > lenp2 .OR. s > lens) THEN
! end of either pattern2 or string
       EXIT ! .FALSE.
     ENDIF
   ENDDO

END FUNCTION match_wild

!PROGRAM test_match
!   IMPLICIT NONE
!   INTEGER, PARAMETER :: np=20, ns=5
!   CHARACTER :: pattern(np)*8, string(ns)*12
!   INTEGER :: s, p
!   LOGICAL :: match_wild
!   EXTERNAL match_wild
!
!   string = (/ 'a.f90       ', 'a1.f90      ', 'a12.f90     ', &
!               'a.f         ', 'logfile.ftp ' /)
!   pattern = (/ 'a*.f90  ', 'a?*.f90 ', 'a*?.f90 ', 'a?*?.f90', &
!                'a*.f90  ', 'a***.f  ', 'a*?*?*? ', 'a**b**c*', &
!                'a*?*b???', '???     ', '*       ', '?       ', &
!                '        ', '**?**   ', '?*      ', '*?*     ', &
!                '*.f90   ', '*?      ', 'a*??.f90', '?????   ' /)
!
!   WRITE(*,'(t17, 100a9)') string
!   DO p = 1,np
!     WRITE(*, '(a, 100L9)') pattern(p), &
!       (match_wild(pattern(p), string(s)), s=1,ns)
!   ENDDO

!END PROGRAM test_match
