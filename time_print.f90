SUBROUTINE time_print
!
! TIME_PRINT Print out the current date and time.
! Notes:
! The standard Fortran 90 routine DATE_AND_TIME is used to get
! the current date and time strings.
!
IMPLICIT NONE

! Data dictionary: declare parameter

INTEGER, PARAMETER :: output = 6

! Data dictionary: declare local scalars.

CHARACTER ( LEN = 8 ) :: datstr
CHARACTER ( LEN = 10 ) :: timstr

! . Get the current date and time.

CALL DATE_AND_TIME ( datstr, timstr )

! . Write out the date and time.

WRITE ( output, "(/A)" ) " Date = " // datstr(7:8) // "/" // &
datstr(5:6) // "/" // &
datstr(1:4)

WRITE ( output, "(A)" ) " Time = " // timstr(1:2) // ":"  // &
timstr(3:4) // ":" // &
timstr(5:10)

WRITE ( output, *)

END SUBROUTINE time_print
