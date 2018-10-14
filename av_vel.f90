SUBROUTINE av_vel ( ntotal, mass, niac, pair_i, pair_j , w, vx, rho, av )

! Subroutine to calculate the average velocity to correct velocity for preventing penetration (monaghan, 1992)

! ntotal : Number of particles                          [in]
! mass   : Particle masses                              [in]
! niac   : Number of interaction pairs                  [in]
! pair_i : List of first partner of interaction pair    [in]
! pair_j : List of second partner of interaction pair   [in]
! w      : Kernel for all interaction pairs             [in]
! vx     : Velocity of each particle                    [in]
! rho    : Density of each particle                     [in]
! av     : Average velocityof each particle            [out]

IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameters types
INTEGER, INTENT(IN) ::  ntotal
INTEGER, INTENT(IN) ::  niac
INTEGER, INTENT(IN) ::  pair_i(max_interaction)
INTEGER, INTENT(IN) ::  pair_j(max_interaction)

REAL(KIND=8), INTENT(IN) :: mass(maxn)
REAL(KIND=8), INTENT(IN) :: w(max_interaction)
REAL(KIND=8), INTENT(IN) :: vx(dim,maxn)
REAL(KIND=8), INTENT(IN) :: rho(maxn)
REAL(KIND=8), INTENT(OUT) :: av(dim, maxn)

! Data dictionary: declare local variable types
INTEGER :: i,j,k,d
REAL (KIND=8) :: vcc, dvx(dim), epsilon

! epsilon --- a small constants chosen by experence, may lead to instability.
! for example, for the 1 dimensional shock tube problem, the E <= 0.3

epsilon = 0.3

! Do Loop 1
DO i = 1, ntotal
   
   ! Do Inner-Loop 1.1
     DO d = 1, dim
        av(d,i) = 0.
     END DO                                       ! End Do Inner-Loop 1.1

END DO                                            ! End Do Loop 1

! Do Loop 2
DO k=1,niac
   i = pair_i(k)
   j = pair_j(k)
   
   ! Do Inner-Loop 2.1
     DO d=1,dim
        dvx(d) = vx(d,i) - vx(d,j)
        av(d, i) = av(d,i) - 2*mass(j)*dvx(d)/(rho(i)+rho(j))*w(k)
        av(d, j) = av(d,j) + 2*mass(i)*dvx(d)/(rho(i)+rho(j))*w(k)
     END DO                                       ! End Do Inner-Loop 2.1

END DO                                            ! End Do Loop 2

! Do Loop 3
DO i = 1, ntotal
   
   ! Do Inner-Loop 3.1
     DO d = 1, dim
        av(d,i) = epsilon * av(d,i)
     END DO                                       ! End Do Inner-Loop 3.1

END DO                                            ! End Do Loop 3

END SUBROUTINE av_vel
