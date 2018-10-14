SUBROUTINE viscosity ( ntotal, itype, x, rho, eta )
!
! Subroutine to define the fluid particle viscosity
!
! ntotal : Number of particles                   [in]
! itype  :  Type of particle                     [in]
! x      : Coordinates of all particles          [in]
! rho    : Density                               [in]
! eta    : Dynamic viscosity                    [out]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

INTEGER, INTENT(IN) :: ntotal
INTEGER, INTENT(IN) :: itype(maxn)

REAL(KIND=8), INTENT(IN) :: x(dim,maxn)
REAL(KIND=8), INTENT(IN) :: rho(maxn)
REAL(KIND=8), INTENT(OUT) :: eta(maxn)

! Data dictionary: declare local variable types

INTEGER :: i

! Do Loop 1
  DO i=1,ntotal

!    If Inner-Loop 1.1
     IF (abs(itype(i)).EQ.1) THEN

        eta(i)=0.

!    Continue Inner-Loop 1.1
     ELSE IF (abs(itype(i)).EQ.2) THEN

        eta(i)=1.0e-3

     END IF                             ! End If Inner-Loop 1.1

  END DO                                ! End Do Loop 1

END SUBROUTINE viscosity
