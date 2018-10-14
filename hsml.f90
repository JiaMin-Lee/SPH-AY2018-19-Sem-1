SUBROUTINE h_upgrade ( dt, ntotal, mass, vx, rho, niac, pair_i, pair_j, dwdx, hsml )
!
! Subroutine to evolve smoothing length
! dt     : time step                                                            [in]
! ntotal : Number of particles                                                  [in]
! mass   : Particle masses                                                      [in]
! vx     : Velocities of all particles                                          [in]
! rho    : Density                                                              [in]
! niac   : Number of interaction pairs                                          [in]
! pair_i : List of first partner of interaction pair                            [in]
! pair_j : List of second partner of interaction pair                           [in]
! dwdx   : Derivative of kernel with respect to x, y and z                      [in]
! hsml   : Smoothing Length                                                 [in/out]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

INTEGER, INTENT(IN) ::  ntotal
INTEGER, INTENT(IN) ::  niac
INTEGER, INTENT(IN) ::  pair_i(max_interaction)
INTEGER, INTENT(IN) ::  pair_j(max_interaction)

REAL(KIND=8), INTENT(IN) :: mass(maxn)
REAL(KIND=8), INTENT(IN) :: vx(dim, maxn)
REAL(KIND=8), INTENT(IN) :: rho(maxn)
REAL(KIND=8), INTENT(IN) :: dwdx(dim, max_interaction)
REAL(KIND=8), INTENT(INOUT) :: hsml(maxn)
REAL(KIND=8), INTENT(IN) :: dt

! Data dictionary: declare local variable types

INTEGER :: i,j,k,d
REAL(KIND=8) :: fac, dvx(dim), hvcc, vcc(maxn), dhsml(maxn)

! If Loop 1 (Block 1)
  IF (sle.EQ.0 ) THEN  

!    --- Keep smoothing length unchanged.
     
     RETURN                                     ! Return to calling program

! Continue Loop 1 (Block 2) 
  ELSE IF (sle.EQ.2) THEN                     

!         --- dh/dt = (-1/dim)*(h/rho)*(drho/dt).
          
!         Do Inner-Loop 1.1
          DO i=1,ntotal
             vcc(i) = 0.e0
          END DO                                ! End Do Inner-Loop 1.1

!         Do Inner-Loop 1.2
          DO k=1,niac
             i = pair_i(k)
             j = pair_j(k)

!            Do Inner-Inner-Loop 1.2.1
             DO d=1,dim
                dvx(d) = vx(d,j) - vx(d,i)
             END DO                             ! End Do Inner-Inner-Loop 1.2.1

             hvcc = dvx(1)*dwdx(1,k)

!            Do Inner-Inner-Loop 1.2.2
             DO d=2,dim
                hvcc = hvcc + dvx(d)*dwdx(d,k)
             END DO                             ! End Do Inner-Inner Loop 1.2.2

             vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
             vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)

          END DO                                ! End Do Inner-Loop 1.2

!         Do Inner-Loop 1.3
          DO i = 1, ntotal
             dhsml(i) = (hsml(i)/dim)*vcc(i)
             hsml(i) = hsml(i) + dt*dhsml(i)

!            If Inner-Inner-Loop 1.3.1
             IF (hsml(i).LE.0) THEN
                hsml(i) = hsml(i) - dt*dhsml(i)
             END IF                             ! End If Inner-Inner-Loop 1.3.1

          END DO                                ! End Do Inner-Loop 1.3

! Continue Loop 1 (Block 3)
  ELSE IF(sle.EQ.1) THEN
         fac = 2.0

!        Do Inner-Loop 1.4
         DO i = 1, ntotal
            hsml(i) = fac * (mass(i)/rho(i))**(1./dim)
         END DO                                 ! End Do Inner-Loop 1.4

  END IF                                        ! End If Loop 1

END SUBROUTINE h_upgrade
