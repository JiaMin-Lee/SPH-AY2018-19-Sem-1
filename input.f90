! Subroutine input

SUBROUTINE input ( x, vx, mass, rho, p, u, itype, hsml, ntotal )
!
! Subroutine for loading or generating initial particle information
! x      : coordinates of particles                      [out]
! vx     : velocities of particles                       [out]
! mass   : mass of particles                             [out]
! rho    : dnesities of particles                        [out]
! p      : pressure of particles                         [out]
! u      : internal energy of particles                  [out]
! itype  : types of particles                            [out]
! hsml   : smoothing lengths of particles                [out]
! ntotal : total particle number                         [out]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data declaration: declare calling parameter types

INTEGER, INTENT(OUT) :: itype(maxn)
INTEGER, INTENT(OUT) :: ntotal

REAL(KIND=8), INTENT(OUT) :: x(dim, maxn)
REAL(KIND=8), INTENT(OUT) :: vx(dim, maxn)
REAL(KIND=8), INTENT(OUT) :: mass(maxn)
REAL(KIND=8), INTENT(OUT) :: p(maxn)
REAL(KIND=8), INTENT(OUT) :: u(maxn)
REAL(KIND=8), INTENT(OUT) :: hsml(maxn)
REAL(KIND=8), INTENT(OUT) :: rho(maxn)

! Data declaration: declare local variable types

INTEGER :: i, d, im

! Load initial particle information from external disk file

! If Loop 1 (Block 1)
  IF (config_input) THEN
    
     OPEN(1,file="../data/f_xv.dat")
     OPEN(2,file="../data/f_state.dat")
     OPEN(3,file="../data/f_other.dat")

     WRITE(*,*)'  **************************************************'
     WRITE(*,*)'      Loading initial particle configuration... '
     READ (1,*) ntotal
     WRITE(*,*)'      Total number of particles    ', ntotal
     WRITE(*,*)'  **************************************************'

!    Do Inner-Loop 1.1
     DO i = 1, ntotal

           read(1,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)
           read(2,*)im, mass(i), rho(i), p(i), u(i)
           read(3,*)im, itype(i), hsml(i)
        
     END DO                                     ! End Do Inner-Loop 1.1

! Continue If Loop 1 (Block 2)
  ELSE

     OPEN(1,file="../data/ini_xv.dat")
     OPEN(2,file="../data/ini_state.dat")
     OPEN(3,file="../data/ini_other.dat")

! Continue If Loop 1 (Block 3)
  IF (shocktube) call shock_tube (x, vx, mass, rho, p, u, itype, hsml, ntotal)
  IF (shearcavity) call shear_cavity (x, vx, mass, rho, p, u, itype, hsml, ntotal)

!    Do Inner-Loop 1.2
     DO i = 1, ntotal

        WRITE(1,1001) i, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)
        WRITE(2,1002) i, mass(i), rho(i), p(i), u(i)
        WRITE(3,1003) i, itype(i), hsml(i)

     END DO                                     ! End Do Inner-Loop 1.2

1001 FORMAT(1x, I6, 6(2x, E15.8))
1002 FORMAT(1x, I6, 7(2x, E15.8))
1003 FORMAT(1x, I6, 2x, I4, 2x, E15.8)

WRITE(*,*)'  **************************************************'
write(*,*)'      Initial particle configuration generated '
write(*,*)'      Total number of particles ', ntotal
WRITE(*,*)'  **************************************************'

  END IF                                        ! End If Loop 1

CLOSE(1)
CLOSE(2)
CLOSE(3)

END SUBROUTINE input

!-------------------------------------------------------------------------------------------------------------------

! Subroutine shock_tube

SUBROUTINE shock_tube ( x, vx, mass, rho, p, u, itype, hsml, ntotal )
!
! This subroutine is used to generate initial data for the
! 1d noh shock tube problem
! x      : coordinates of particles                     [out]
! vx     : velocities of particles                      [out]
! mass   : mass of particles                            [out]
! rho    : densities of particles                       [out]
! p      : pressure of particles                        [out]
! u      : internal energy of particles                 [out]
! itype  : types of particles                           [out]
!          =1 ideal gas
! hsml   : smoothing lengths of particles               [out]
! ntotal : total particle number                        [out]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

INTEGER, INTENT(OUT) :: itype(maxn)
INTEGER, INTENT(OUT) :: ntotal

REAL(KIND=8), INTENT(OUT) :: x(dim, maxn)
REAL(KIND=8), INTENT(OUT) :: vx(dim, maxn)
REAL(KIND=8), INTENT(OUT) :: mass(maxn)
REAL(KIND=8), INTENT(OUT) :: rho(maxn)
REAL(KIND=8), INTENT(OUT) :: p(maxn)
REAL(KIND=8), INTENT(OUT) :: u(maxn)
REAL(KIND=8), INTENT(OUT) :: hsml(maxn)

! Data dictionary: declare local variables

INTEGER :: i, d

REAL(KIND=8) :: space_x

ntotal=400
space_x=0.6/80.

! Do Loop 1
  DO i=1,ntotal

     mass(i)=0.75/400.
     hsml(i)=0.015
     itype(i)=1

!    Do Inner-Loop 1.1
     DO d = 1, dim

        x(d,i) = 0.
        vx(d,i) = 0.

     END DO                             ! End Do Inner-Loop 1.1

  END DO                                ! End Do Loop 1

! Do Loop 2
  DO i=1,320

     x(1,i)=-0.6+space_x/4.*(i-1)

  END DO                                ! End Do Loop 2

! Do Loop 3
  DO i=320+1,ntotal

     x(1,i)=0.+space_x*(i-320)

  END DO                                ! End Do Loop 3

! Do Loop 4
  DO i=1,ntotal

!    If Inner-Loop 4.1
     IF (x(1,i).LE.1.E-8) THEN

        u(i)=2.5
        rho(i)=1.
        p(i)=1.

     END IF                             ! End If Inner-Loop 4.1

!    If Inner-Loop 4.2
     IF (x(1,i).GT.1.E-8) THEN

        u(i)=1.795
        rho(i)=0.25
        p(i)=0.1795

     END IF                             ! End If Inner-Loop 4.2

  END DO                                ! End Do Loop 4

END SUBROUTINE shock_tube

!------------------------------------------------------------------------------------------------------------------

! Subroutine shear cavity

SUBROUTINE shear_cavity ( x, vx, mass, rho, p, u, itype, hsml, ntotal )
!
! This subroutine is used to generate initial data for the
! 2 d shear driven cavity probem with Re = 1
! x      : coordinates of particles                             [out]
! vx     : velocities of particles                              [out]
! mass   : mass of particles                                    [out]
! rho    : dnesities of particles                               [out]
! p      : pressure of particles                                [out]
! u      : internal energy of particles                         [out]
! itype  : types of particles                                   [out]
!          =2 water
! h      : smoothing lengths of particles                       [out]
! ntotal : total particle number                                [out]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

INTEGER, INTENT(OUT) :: itype(maxn)
INTEGER, INTENT(OUT) :: ntotal

REAL(KIND=8), INTENT(OUT) :: x(dim, maxn)
REAL(KIND=8), INTENT(OUT) :: vx(dim,maxn)
REAL(KIND=8), INTENT(OUT) :: mass(maxn)
REAL(KIND=8), INTENT(OUT) :: rho(maxn)
REAL(KIND=8), INTENT(OUT) :: p(maxn)
REAL(KIND=8), INTENT(OUT) :: u(maxn)
REAL(KIND=8), INTENT(OUT) :: hsml(maxn)

! Data dictionary: declare local variable types

INTEGER :: i, j, d, m, n, mp, np, k

REAL(KIND=8) :: x1, y1, dx, dy

! Giving mass and smoothing length as well as other data.

  m = 41
  n = 41
  mp = m-1
  np = n-1
  ntotal = mp * np
  x1 = 1.E-3
  y1 = 1.E-3
  dx = x1/mp
  dy = y1/np

! Do Loop 1
  DO i = 1, mp

!    Do Inner-Loop 1.1
     DO j = 1, np

        k = j + (i-1)*np
        x(1,k) = (i-1)*dx + dx/2.
        x(2,k) = (j-1)*dy + dy/2.

     END DO                             ! End Do Inner-Loop 1.1

  END DO                                ! End Do Loop 1

! Do Loop 2
  DO i = 1, mp*np

     vx(1, i) = 0.
     vx(2, i) = 0.
     rho (i) = 1000.
     mass(i) = dx*dy*rho(i)
     p(i)= 0.
     u(i)=357.1
     itype(i) = 2
     hsml(i) = dx

  END DO                                ! End Do Loop 2

END SUBROUTINE shear_cavity
