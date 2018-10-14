SUBROUTINE ext_force ( ntotal, mass, x, niac, pair_i, pair_j, itype, hsml, dvxdt )
!
! Subroutine to calculate the external forces, e.g. gravitational forces.
! The forces from the interactions with boundary virtual particles are also calculated here as external forces.
!
! here as the external force.
! ntotal : Number of particles                                  [in]
! mass   : Particle masses                                      [in]
! x      : Coordinates of all particles                         [in]
! pair_i : List of first partner of interaction pair            [in]
! pair_j : List of second partner of interaction pair           [in]
! itype  : type of particles                                    [in]
! hsml   : Smoothing Length                                     [in]
! dvxdt  : Acceleration with respect to x, y and z             [out]
!
IMPLICIT NONE
INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameters types

INTEGER, INTENT(IN) ::  ntotal
INTEGER, INTENT(IN) ::  itype(maxn)
INTEGER, INTENT(IN) ::  niac
INTEGER, INTENT(IN) ::  pair_i(max_interaction)
INTEGER, INTENT(IN) ::  pair_j(max_interaction)

REAL(KIND=8), INTENT(IN) :: mass(maxn)
REAL(KIND=8), INTENT(IN) :: x(dim,maxn)
REAL(KIND=8), INTENT(IN) :: hsml(maxn)
REAL(KIND=8), INTENT(OUT) :: dvxdt (dim,maxn)

! Data Dictionary: declare local variable types

INTEGER i, j, k, d

REAL(KIND=8) :: dx(dim), rr, f, rrO, dd, p1, p2

! Do Loop 1
  DO i = 1, ntotal

!    Do Inner-Loop 1.1
     DO d = 1, dim
        dvxdt(d, i) = 0.
     END DO                                     ! End Do Inner-Loop 1.1

  END DO                                        ! End Do Loop 1

! Consider self-gravity or not ?
! If Loop 2
  IF (self_gravity) THEN

!    Do Inner-Loop 2.1
     DO i = 1, ntotal
        dvxdt(dim, i) = -9.8
     END DO                                     ! End Do Inner-Loop 2.1

  END IF                                        ! End If Loop 2

! Boundary particle force and penalty anti-penetration force.
rrO = 1.25e-5
dd = 1.e-2
p1 = 12
p2 = 4

! Do Loop 3
  DO k=1,niac
     i = pair_i(k)
     j = pair_j(k)

!    If Inner-Loop 3.1
     IF(itype(i).GT.0.AND.itype(j).LT.0) THEN
        rr = 0.

!       Do Inner-Inner-Loop 3.1.1
        DO d=1,dim
           dx(d) = x(d,i) - x(d,j)
           rr = rr + dx(d)*dx(d)
        END DO                                  ! End Do Inner-Inner-Loop 3.1.1

        rr = sqrt(rr)

!       If Inner-Inner-Loop 3.1.2
        IF(rr.LT.rrO) THEN
          f = ((rrO/rr)**p1-(rrO/rr)**p2)/rr**2

!         Do Inner-Inner-Inner-Loop 3.1.2.1
          DO d = 1, dim
             dvxdt(d, i) = dvxdt(d, i) + dd*dx(d)*f
          END DO                                ! End Do Inner-Inner-Inner-Loop 3.1.2.1

        END IF                                  ! End If Inner-Inner Loop 3.1.2

     END IF                                     ! End If Inner-Loop 3.1

  END DO                                        ! End Do Loop 3

END SUBROUTINE ext_force
