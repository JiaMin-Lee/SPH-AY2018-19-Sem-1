PROGRAM SPH
!
! This is a three dimensional SPH code, the followings are the
! basic parameters needed in this codeor calculated by this code
!
! mass   : mass of particles                             [in]
! ntotal : total particle number ues                     [in]
! dt     : Time step used in the time integration        [in]
! itype  : types of particles                            [in]
! x      : coordinates of particles                  [in/out]
! vx     : velocities of particles                   [in/out]
! rho    : dnesities of particles                    [in/out]
! p      : pressure of particles                     [in/out]
! u      : internal energy of particles              [in/out]
! hsml   : smoothing lengths of particles            [in/out]
! c      : sound velocity of particles                  [out]
! s      : entropy of particles                         [out]
! e      : total energy of particles                    [out]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

INTEGER :: ntotal
INTEGER :: itype(maxn)
!INTEGER, INTENT(OUT) :: maxtimestep                  !included

REAL(KIND=8) :: x(dim,maxn)
REAL(KIND=8) :: vx(dim, maxn)
REAL(KIND=8) :: mass(maxn)
REAL(KIND=8) :: rho(maxn)
REAL(KIND=8) :: p(maxn)
REAL(KIND=8) :: u(maxn)
REAL(KIND=8) :: c(maxn)
REAL(KIND=8) :: s(maxn)
REAL(KIND=8) :: e(maxn)
REAL(KIND=8) :: hsml(maxn)
REAL(KIND=8) :: dt

! Data dictionary: declare local variable types

INTEGER :: maxtimestep, d, m, i, yesorno
INTEGER :: numx, numy                           !included
!INTEGER ::  d, m, i, yesorno                   !included
REAL(KIND=8) :: s1, s2


CALL time_print
CALL time_elapsed(s1)


IF (shocktube) dt = 0.005
IF (shearcavity) dt = 5.e-5

CALL input ( x, vx, mass, rho, p, u, itype, hsml, ntotal )

1 WRITE(*,*)' ***************************************************'
  WRITE(*,*)' Please input the maximal time steps '
  WRITE(*,*)' ***************************************************'
  READ(*,*) maxtimestep
  WRITE(*,*)' For interpolation into mesh, please select number of grid in x-direction '
  READ(*,*) numx
  WRITE(*,*)' For interpolation into mesh, please select number of grid in y-direction '
  READ(*,*) numy

CALL time_integration ( x, vx, mass, rho, p, u, c, s, e, itype, hsml, ntotal, maxtimestep, numx, numy, dt )
CALL output ( x, vx, mass, rho, p, u, c, itype, hsml, ntotal, maxtimestep, numx, numy,  dt )       !included dt,num  and maxtimestep


WRITE(*,*)' ***************************************************'
WRITE(*,*) 'Are you going to run more time steps ? (0=No, 1=yes)'
READ (*,*) yesorno

IF (yesorno.NE.0) GO TO 1


CALL time_print
CALL time_elapsed(s2)

WRITE (*,*)'         Elapsed CPU time = ', s2-s1

END PROGRAM SPH
