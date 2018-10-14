SUBROUTINE p_gas ( rho, u, p, c )
!
! Gamma law EOS: subroutine to calculate the pressure and sound
! rho : Density                  [in]
! u   : Internal energy          [in]
! p   : Pressure                [out]
! c   : sound velocity          [out]
!
IMPLICIT NONE

! Data dictionary: declare calling parameter types

REAL(KIND=8), INTENT(IN) :: rho
REAL(KIND=8), INTENT(IN) :: u
REAL(KIND=8), INTENT(OUT) :: p
REAL(KIND=8), INTENT(OUT) :: c
!REAL(KIND=8), INTENT(OUT) :: gamma                      !included

! Data dictionary: declare local variable types

REAL(KIND=8) :: gamma

! For air (idea gas)

  gamma=1.4
  p = (gamma-1) * rho * u
  c = sqrt((gamma-1) * u)

END SUBROUTINE p_gas

!-------------------------------------------------------------------------------------------------------------------

SUBROUTINE p_art_water ( rho, p, c, gamma)
!
! Artificial equation of state for the artificial compressibility
! rho : Density                  [in]
! u   : Internal energy          [in]
! p   : Pressure                [out]
! c   : sound velocity          [out]
! Equation of state for artificial compressibility
!
IMPLICIT NONE

! Data dicionary: declare calling parameter types

REAL(KIND=8), INTENT(IN) :: rho
REAL(KIND=8), INTENT(OUT) :: p 
REAL(KIND=8), INTENT(OUT) :: c
!REAL(KIND=8), INTENT(OUT) :: gamma

! Data dicionary: declare local variable types

REAL(KIND=8) :: rhoO, u
REAL(KIND=8) :: gamma

! Artificial EOS, Form 1 (Monaghan, 1994)
! gamma=7.
! rho0=1000.
! b = 1.013e5
! p = b*{(rho/rhoO)**gamma-1)
! c = 1480.

! Artificial EOS, Form 2 (Morris, 1997)
  c = 0.01
  p = c**2 * rho

END SUBROUTINE p_art_water
