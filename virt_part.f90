SUBROUTINE virt_part ( itimestep, ntotal, nvirt, hsml, mass, x, vx, rho, u, p, itype )
!
! Subroutine to determine the information of virtual particles
! Here only the Monaghan type virtual particles for the 2D shear
! cavity driven problem are generated.
!
! itimestep : Current time step                             [in]
! ntotal    : Number of particles                           [in]
! nvirt     : Number of virtual particles                  [out]
! hsml      : Smoothing Length                          [in/out]
! mass      : Particle masses                           [in/out]
! x         : Coordinates of all particles              [in/out]
! vx        : Velocities of all particles               [in/out]
! rho       : Density                                   [in/out]
! u         : internal energy                           [in/out]
! itype     : type of particles                         [in/out]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

INTEGER, INTENT(IN) :: itimestep
INTEGER, INTENT(IN) :: ntotal
INTEGER, INTENT(INOUT) :: itype(maxn)
INTEGER, INTENT(OUT) :: nvirt

REAL(KIND=8), INTENT(INOUT) :: hsml(maxn) 
REAL(KIND=8), INTENT(INOUT) :: mass (maxn) 
REAL(KIND=8), INTENT(INOUT) :: x(dim,maxn) 
REAL(KIND=8), INTENT(INOUT) :: vx(dim,maxn)
REAL(KIND=8), INTENT(INOUT) :: rho (maxn)
REAL(KIND=8), INTENT(INOUT) :: u(maxn)
REAL(KIND=8), INTENT(INOUT) :: p(maxn)

! Data dictionary: declare local variable types

INTEGER :: i, j, d, im, mp

REAL(KIND=8) :: x1, dx, v_inf

! If Loop 1
  IF (vp_input) THEN

     OPEN(1,file="../data/xv_vp.dat")
     OPEN(2,file="../data/state_vp.dat")
     OPEN(3,file="../data/other_vp.dat")

     READ(1,*) nvirt

!    DO Inner-Loop 1.1
     DO j = 1, nvirt

        i = ntotal + j
        
        READ(1,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)
        READ(2,*)im, mass(i), rho(i), p(i), u(i)
        READ(3,*)im, itype(i), hsml(i)

     END DO                             ! End Do Inner-Loop 1.1

     CLOSE(1)
     CLOSE(2)
     CLOSE(3)

! Continue If Loop 1
  ELSE

     nvirt = 0
     mp = 40 
     x1 = 1.0e-3
     dx = x1 / mp
     v_inf = 1.e-3

! Monaghan type virtual particle on the Upper side

!    Do Inner-Loop 1.2
     DO i = 1, 2*mp+1

        nvirt = nvirt + 1
        x(1, ntotal + nvirt) = (i-1)*dx/2
        x(2, ntotal + nvirt) = x1
        vx(1, ntotal + nvirt) = v_inf
        vx(2, ntotal + nvirt) = 0.

     END DO                             ! End Do Inner-Loop 1.2

! Monaghan type virtual particle on the Lower side

!    Do Inner-Loop 1.3
     DO i = 1, 2*mp+1

        nvirt = nvirt + 1
        x(1, ntotal + nvirt) = (i-1)*dx/2
        x(2, ntotal + nvirt) = 0.
        vx(1, ntotal + nvirt) = 0.
        vx(2, ntotal + nvirt) = 0.

     END DO                             ! End Do Inner-Loop 1.3

! Monaghan type virtual particle on the Left side

!    Do Inner-Loop 1.4
     DO i = 1, 2*mp-1

        nvirt = nvirt + 1
        x(1, ntotal + nvirt) = 0.
        x(2, ntotal + nvirt) = i*dx/2
        vx(1, ntotal + nvirt) = 0.
        vx(2, ntotal + nvirt) = 0.

     END DO                             ! End Do Inner-Loop 1.4

! Monaghan type virtual particle on the Right side

!    Do Inner-Loop 1.5
     DO i = 1, 2*mp-1

        nvirt = nvirt + 1
        x(1, ntotal + nvirt) = x1
        x(2, ntotal + nvirt) = i*dx/2
        vx(1, ntotal + nvirt) = 0.
        vx(2, ntotal + nvirt) = 0.

     END DO                             ! End Do Inner-Loop 1.5

!    Do Inner-Loop 1.6
     DO i = 1, nvirt

        rho (ntotal + i) = 1000.
        mass(ntotal + i) = rho (ntotal + i) * dx * dx
        p(ntotal + i) = 0.
        u(ntotal + i) = 357.1
        itype(ntotal + i) = -2
        hsml(ntotal + i) = dx

     END DO                             ! End Do Inner-Loop 1.6

  END IF                                ! End If Loop 1

! If Loop 2
  IF (mod (itimestep, save_step) .EQ.0) THEN

     OPEN(1,file="../data/xv_vp.dat")
     OPEN(2,file="../data/state_vp.dat")
     OPEN(3,file="../data/other_vp.dat")

     !WRITE(1,*) nvirt                                  !excluded

!    Do Inner-Loop 2.1
     DO i = ntotal + 1, ntotal + nvirt

        WRITE(1,1001) i, (x(d, i) , d=1,dim), (vx(d,i) , d = 1, dim)
        WRITE(2,1002) i, mass(i), rho(i), p(i), u(i)
        WRITE(3,1003) i, itype(i), hsml(i)

     END DO                             ! End Do Inner-Loop 2.1

    1001 FORMAT(1x, I6, 6(2x, e15.8))
    1002 FORMAT(1x, I6, 7(2x, e15.8))
    1003 FORMAT(1x, I6, 2x, I4, 2x, e15.8)

    CLOSE(1)
    CLOSE(2)
    CLOSE(3)

  END IF                                ! End If Loop 

! If Loop 3
  IF (mod(itimestep,print_step).EQ.0) THEN

!    If Inner-Loop 3.1
     IF (int_stat) THEN

        WRITE (*,*) ' >> Statistics: Virtual boundary particles:'
        WRITE (*,*) ' Number of virtual particles:',NVIRT

     END IF                             ! End If Inner-Loop 3.1

  END IF                                ! End If Loop 3

END SUBROUTINE virt_part
