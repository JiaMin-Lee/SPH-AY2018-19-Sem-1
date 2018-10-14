! Subroutine sum_density
SUBROUTINE sum_density ( ntotal, hsml, mass, niac, pair_i, pair_j , w, itype, rho )
!
! Subroutine to calculate the density with SPH summation algorithm.
! ntotal : Number of particles                                  [in]
! hsml   : Smoothing Length                                     [in]
! mass   : Particle masses                                      [in]
! niac   : Number of interaction pairs                          [in]
! pair_i : List of first partner of interaction pair            [in]
! pair_j : List of second partner of interaction pair           [in]
! w      : Kernel for all interaction pairs                     [in]
! itype  : type of particles                                    [in]
! x      : Coordinates of all particles                         [in]
! rho    : Density                                             [out]

IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types
INTEGER, INTENT(IN) :: ntotal
INTEGER, INTENT(IN) :: niac
INTEGER, INTENT(IN) :: pair_i(max_interaction)
INTEGER, INTENT(IN) :: pair_j(max_interaction)
INTEGER, INTENT(IN) :: itype(maxn)

REAL(KIND=8), INTENT(IN) :: hsml(maxn)
REAL(KIND=8), INTENT(IN) :: mass(maxn)
REAL(KIND=8), INTENT(IN) :: w(max_interaction)
REAL(KIND=8), INTENT(OUT) :: rho (maxn)

! Data dictionary: declare local variable types
INTEGER :: i, j, k, d

REAL(KIND=8) :: selfdens, hv(dim), r, wi(maxn)

! wi(maxn) --- integration of the kernel itself

!  Do Loop 1
   DO d=1,dim
      hv(d) = 0.0D0
   END DO                                         ! End Do Loop 1
 
! Self density of each particle: Wii (Kernel for distance 0) and take contribution of particle itself:

  r=0.

! Firstly, calculate the integration of the kernel over the space

!  Do Loop 2
   DO i=1,ntotal
      ! Call subroutine kernel
        CALL kernel (r, hv, hsml(i), selfdens, hv)
     
      wi(i)=selfdens*mass(i)/rho(i)
   END DO                                         ! End Do Loop 2

!  Do Loop 3
   DO k=1,niac
      i = pair_i(k)
      j = pair_j(k)
      wi(i) = wi(i) + mass(j)/rho(j)*w(k)
      wi(j) = wi(j) + mass(i)/rho(i)*w(k)
   END DO                                         ! End Do Loop 3

! Secondly calculate the rho integration over the space

!  Do Loop 4
   DO i=1,ntotal
      ! Call subroutine kernel
        CALL kernel(r, hv, hsml(i), selfdens, hv)

      rho(i) = selfdens*mass(i)
   END DO                                         ! End Do Loop 4

! Calculate SPH sum for rho:

!  Do Loop 5
   DO k=1,niac
      i = pair_i(k)
      j = pair_j(k)
      rho(i) = rho(i) + mass(j)*w(k)
      rho(j) = rho(j) + mass(i)*w(k)
   END DO                                         ! End Do Loop 5

! Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)

!  If Loop 6
   IF (nor_density) THEN
   ! Do Inner-Loop 6.1
     DO i=1, ntotal
        rho(i)=rho(i)/wi(i)
     END DO                                       ! End Do Inner-Loop 6.1
   END IF                                         ! End If Loop 6

END SUBROUTINE sum_density

!------------------------------------------------------------------------------------------------------------------

! Subroutine con_density
SUBROUTINE con_density ( ntotal, mass, niac, pair_i, pair_j, dwdx, vx, itype, x, rho, drhodt )
!
! Subroutine to calculate 'the density with SPH continuiity approach.
! ntotal : Number of particles                                  [in]
! mass   : Particle masses                                      [in]
! niac   : Number of interaction pairs                          [in]
! pair_i : List of first partner of interaction pair            [in]
! pair_j : List of second partner of interaction pair           [in]
! dwdx   : derivation of Kernel for all interaction pairs       [in]
! vx     : Velocities of all particles                          [in]
! itype  : type of particles                                    [in]
! x      : Coordinates of all particles                         [in]
! rho    : Density                                              [in]
! drhodt : Density change rate of each particle                [out]

IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

INTEGER, INTENT(IN) :: ntotal
INTEGER, INTENT(IN) :: niac
INTEGER, INTENT(IN) :: pair_i(max_interaction)
INTEGER, INTENT(IN) :: pair_j(max_interaction)
INTEGER, INTENT(IN) :: itype(maxn)

REAL(KIND=8), INTENT(IN) :: mass(maxn)
REAL(KIND=8), INTENT(IN) :: dwdx(dim, max_interaction)
REAL(KIND=8), INTENT(IN) :: vx(dim,maxn)
REAL(KIND=8), INTENT(IN) :: x(dim,maxn)
REAL(KIND=8), INTENT(IN) :: rho(maxn)
REAL(KIND=8), INTENT(OUT) :: drhodt(maxn)

! Data dictionary: declare local variable types

INTEGER :: i,j,k,d

REAL(KIND=8) :: vcc, dvx(dim)

! Do Loop 1
  DO i = 1, ntotal
     drhodt(i) = 0.
  END DO                                          ! End Do Loop 1

! Do Loop 2
  DO k=1,niac
     i = pair_i(k)
     j = pair_j(k)

  ! Do Inner-Loop 2.1
    DO d=1,dim
       dvx(d) = vx(d,i) - vx(d,j)
    END DO                                        ! End Do Inner-Loop 2.1

  vcc = dvx(1)*dwdx(1,k)

  ! Do Inner-Loop 2.2
    DO d=2, dim
       vcc = vcc + dvx(d)*dwdx(d,k)
    END DO                                        ! End Do Inner-Loop 2.2

  drhodt(i) = drhodt(i) + mass(j)*vcc
  drhodt(j) = drhodt(j) + mass(i)*vcc

  END DO                                          ! End Do Loop 2

END SUBROUTINE con_density
