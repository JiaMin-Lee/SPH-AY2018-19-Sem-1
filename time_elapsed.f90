SUBROUTINE time_elapsed ( s )
!
! The standard Fortran 90 routine RTC is used to calculate the elapsed CPU
!
USE ifport

IMPLICIT NONE

! Data dictionary: declare parameter
INTEGER, PARAMETER :: output = 6

! Data dictionary: declare calling parameter types
REAL(KIND=8), INTENT(OUT) :: s

s = rtc()

END SUBROUTINE time_elapsed
