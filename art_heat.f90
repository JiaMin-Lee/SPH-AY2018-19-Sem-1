SUBROUTINE art_heat ( ntotal, hsml, mass, x, vx, niac, rho, u, c, pair_i, pair_j, w, dwdx, dedt )
!
! Subroutine to calculate the artificial heatlFulk, 1994, p, a-17)
!
! ntotal : Number of particles                                          [in]
! hsml   : Smoothing Length                                             [in]
! mass   : Particle masses                                              [in]
! x      : Coordinates of all particles                                 [in]
! vx     : Velocities of all particles                                  [in]
! rho    : Density                                                      [in]
! u      : specific internal energy                                     [in]
! c      : Sound veolcity                                               [in]
! niac   : Number of interaction pairs                                  [in]
! pair_i : List of first partner of interaction pair                    [in]
! pair_j : List of second partner of interaction pair                   [in]
! w : Kernel for all interaction pairs                                  [in]
! dwdx   : Derivative of kernel with respect to x, y and z              [in]
! dedt   : produced artificial heat, adding to energy Eq.              [out]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types
INTEGER, INTENT(IN) ::  ntotal
INTEGER, INTENT(IN) ::  niac
INTEGER, INTENT(IN) ::  pair_i(max_interaction)
INTEGER, INTENT(IN) ::  pair_j(max_interaction)

REAL(KIND=8), INTENT(IN) :: hsml(maxn)
REAL(KIND=8), INTENT(IN) :: mass(maxn)
REAL(KIND=8), INTENT(IN) :: x(dim,maxn)
REAL(KIND=8), INTENT(IN) :: vx(dim,maxn)
REAL(KIND=8), INTENT(IN) :: rho(maxn)
REAL(KIND=8), INTENT(IN) :: u(maxn)
REAL(KIND=8), INTENT(IN) :: c(maxn)
REAL(KIND=8), INTENT(IN) :: w(max_interaction)
REAL(KIND=8), INTENT(IN) :: dwdx(dim,max_interaction)
REAL(KIND=8), INTENT(OUT) :: dedt(maxn)

! Data dictionary: declare local variable types
INTEGER :: i,j,k,d

REAL(KIND=8) :: dx, dvx(dim), vr, rr, h, mc, mrho, mhsml, vcc(maxn), hvcc, mui, muj, muij, rdwdx, g1, g2

! Parameter for the artificial heat conduction:
g1=0.1
g2=1.0

! Do Loop 1:
DO i=1,ntotal
   vcc(i) = 0.0D0
   dedt(i) = 0.0D0
END DO                                          ! End Do Loop 1

! Do Loop 2:
DO k=1,niac
   i = pair_i(k)
   j = pair_j(k)
   
   ! Do Inner-Loop 2.1:
   DO d=1,dim
      dvx(d) = vx(d,j) - vx(d,i)
   END DO                                        ! End Do Inner-Loop 2.1
   hvcc = dvx(1)*dwdx(1,k)
   
   ! Do Inner-Loop 2.2:
     DO d=2,dim
        hvcc = hvcc + dvx(d)*dwdx(d,k)
     END DO                                      ! End Do Inner-Loop 2.2
   
vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
END DO                                           ! End Do Loop 2

! Do Loop 3:
DO k=1,niac
   i = pair_i(k)
   j = pair_j(k)
   mhsml= (hsml(i)+hsml(j))/2.
   mrho = (0.5D0)*(rho(i) + rho(j))
   rr = 0.0D0
   rdwdx = 0.0D0

   ! Do Inner-Loop 3.1: 
   DO d=1,dim
      dx = x(d,i) - x(d,j)
      rr = rr + dx*dx
      rdwdx = rdwdx + dx*dwdx(d,k)
   END DO                                        ! End Do Inner-Loop 3.1

   mui=g1*hsml(i)*c(i)+g2*hsml(i)**2*(abs(vcc(i))-vcc(i))
   muj=g1*hsml(j)*c(j)+g2*hsml(j)**2*(abs(vcc(j))-vcc(j))
   muij= 0.5*(mui+muj)
   h = muij/(mrho*(rr+0.01*mhsml**2))*rdwdx
   dedt(i) = dedt(i) + mass(j)*h*(u(i)-u(j))
   dedt(j) = dedt(j) + mass(i)*h*(u(j)-u(i))
END DO                                           ! End Do Loop 3

! Do Loop 4:
DO i=1,ntotal
   dedt(i) = 2.0D0*dedt(i)
END DO                                           ! End Do Loop 4

END SUBROUTINE art_heat
