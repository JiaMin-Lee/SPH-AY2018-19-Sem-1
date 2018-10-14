SUBROUTINE art_visc ( ntotal, hsml, mass, x, vx, niac, rho, c, pair_i, pair_j, w, dwdx, dvxdt, dedt )
!
! Subroutine to calculate the artificial viscosity (Monaghan, 1992)
!
! ntotal : Number of particles (including virtual particles)    [in]
! hsml   : Smoothing Length                                     [in]
! mass   : Particle masses                                      [in]
! x      : Coordinates of all particles                         [in]
! vx     : Velocities of all particles                          [in]
! niac   : Number of interaction pairs                          [in]
! rho    : Density                                              [in]
! c      : Temperature                                          [in]
! pair_i : List of first partner of interaction pair            [in]
! pair_j : List of second partner of interaction pair           [in]
! w      : Kernel for all interaction pairs                     [in]
! dwdx   : Derivative of kernel with respect to x, y and z      [in]
! dvxdt  : Acceleration with respect to x, y and z             [out]
! dedt   : Change of specific internal energy                  [out]

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
REAL(KIND=8), INTENT(IN) :: c(maxn)
REAL(KIND=8), INTENT(IN) :: w(max_interaction)
REAL(KIND=8), INTENT(IN) :: dwdx(dim,max_interaction)
REAL(KIND=8), INTENT(OUT) :: dvxdt(dim,maxn)
REAL(KIND=8), INTENT(OUT) :: dedt(maxn)

! Data dictionary: declare local variable types
INTEGER :: i,j,k,d

REAL(KIND=8) :: dx, dvx(dim), piv, muv, vr, rr, h, mc, mrho, mhsml

! Data dictionary: declare constants
! Parameter for the artificial viscosity:
! Shear viscosity
REAL(KIND=8), PARAMETER :: alpha = 1.0D0

! Bulk viscosity
REAL(KIND=8), PARAMETER :: beta = 1.0D0

! Parameter to avoid singularities
REAL(KIND=8), PARAMETER :: etq = 0.1D0

! Do Loop 1
DO i=1,ntotal
   
   ! Do Inner-Loop 1.1
   DO d=1,dim
      dvxdt(d,i) = 0.0D0
   END DO                                              ! End Inner-Loop 1.1

   dedt(i) = 0.0D0
END DO                                                 ! End Do Loop 1

! Calculate SPH sum for artificial viscosity
! Do Loop 2
DO k=1,niac
   i = pair_i(k)
   j = pair_j(k)
   mhsml= (hsml(i)+hsml(j))/2
   vr = 0.0D0
   rr = 0.0D0

   ! Do Inner-Loop 2.1
   DO d=1,dim
      dvx(d) = vx(d,i) - vx(d,j)
      dx = x(d,i) - x(d,j)
      vr = vr + dvx(d)*dx
      rr = rr + dx*dx
   END DO                                               ! End Do Loop 2.1

! Artificial viscous force only if v_ij * r_ij < 0

! If Inner-Loop 2.2 
  IF (vr.LT.0.0D0) THEN
      ! Calculate muv_ij = hsml v_ij * r_ij / ( r_ijA2 + hsml*2 etq*2 )
        muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq)
      ! Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij*2) / rho_ij
        mc = 0.5D0*(c(i) + c(j))
        mrho = 0.5D0*(rho(i) + rho(j))
        piv = (beta*muv - alpha*mc)*muv/mrho

      ! Calculate SPH sum for artificial viscous force
      ! Do Inner-Inner-Loop 2.2.1  
         DO d=1,dim
            h = -piv*dwdx(d,k)
            dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
            dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
            dedt(i) = dedt(i) - mass(j)*dvx(d)*h
            dedt(j) = dedt(j) - mass(i)*dvx(d)*h
         END DO                                          ! End Do Inner-Inner-Loop 2.2.1
  END IF                                                 ! End If Inner-Loop 2.2
END DO                                                   ! End Do Loop 2

! Change of specific internal energy:
! Do Loop 3
DO i=1,ntotal
   dedt(i) = 0.5D0*dedt(i)
END DO                                                   ! End Do Loop 3

END SUBROUTINE art_visc
