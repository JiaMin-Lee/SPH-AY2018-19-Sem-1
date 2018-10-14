SUBROUTINE time_integration ( x, vx, mass, rho, p, u, c, s, e, itype, hsml, ntotal, maxtimestep, dt )
!
! x           : coordinates of particles [input/output]
! vx          : velocities of particles [input/output]
! mass        : mass of particles [input]
! rho         : densities of particles [input/output]
! p           : pressure of particles [input/output]
! u           : internal energy of particles [input/output]
! c           : sound velocity of particles [output]
! s           : entropy of particles, not used here [output]
! e           : total energy of particles [output]
! itype       : types of particles [input]
!               =1 ideal gas
!               =2 water
!               =3 tnt
! hsml        : smoothing lengths of particles [input/output]
! ntotal      : total particle number [input]
! maxtimestep : maximum timesteps [input]
! dt          : timestep [input]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

INTEGER, INTENT(IN) :: itype(maxn)
INTEGER, INTENT(IN) :: ntotal
INTEGER, INTENT(IN) :: maxtimestep

REAL(KIND=8),INTENT(IN) :: mass(maxn)
REAL(KIND=8),INTENT(IN) :: dt
REAL(KIND=8),INTENT(INOUT) :: x(dim, maxn)
REAL(KIND=8),INTENT(INOUT) :: vx(dim, maxn)
REAL(KIND=8),INTENT(INOUT) :: rho(maxn)
REAL(KIND=8),INTENT(INOUT) :: p(maxn)
REAL(KIND=8),INTENT(INOUT) :: u(maxn)
REAL(KIND=8),INTENT(INOUT) :: hsml(maxn)
REAL(KIND=8),INTENT(OUT) :: c(maxn)
REAL(KIND=8),INTENT(OUT) :: s(maxn)
REAL(KIND=8),INTENT(OUT) :: e(maxn)

! Data dictionary: declare local variable types

INTEGER :: i, j, k, itimestep, d, current_ts, nstart

REAL(KIND=8) :: x_min(dim, maxn), v_min(dim, maxn), u_min(maxn),rho_min(maxn), dx(dim,maxn), dvx(dim, maxn), du(maxn)
REAL(KIND=8) :: drho(maxn), av(dim, maxn), ds(maxn), t(maxn), tdsdt(maxn)
REAL(KIND=8) :: time, temp_rho, temp_u

! Do Loop 1
  DO i = 1, ntotal

!    Do Inner-Loop 1.1
     DO d = 1, dim

        av(d, i) = 0.

     END DO                                     ! End Do Inner-Loop 1.1

  END DO                                        ! End Do Loop 1

! Do Loop 2
  DO itimestep = nstart+1, nstart+maxtimestep

     current_ts = current_ts + 1

!    If Inner-Loop 2.1
     IF (mod(itimestep,print_step).EQ.0) THEN

        WRITE(*,*)' ___________________________________________________'
        WRITE(*,*)'     current number of time step =', itimestep,' current time=', real(time+dt)
        WRITE(*,*)' ___________________________________________________'

     END IF                                     ! End If Inner-Loop 2.1

! If not first time step, then update thermal energy, density and velocity half a time step

!    If Inner-Loop 2.2
     IF (itimestep .NE. 1) THEN

!       Do Inner-Inner-Loop 2.2.1
        DO i = 1, ntotal

           u_min(i) = u(i)
           temp_u = 0.
 
           IF (dim.EQ.1) temp_u=-nsym*p(i)*vx(1,i)/x(1,i)/rho(i)

           u(i) = u(i) + (dt/2.)* (du(i)+temp_u)

           IF (u(i).LT.0) u(i) = 0.

!          If Inner-Inner-Inner-Loop 2.2.1.1
           IF (.NOT.summation_density) THEN

              rho_min(i) = rho(i)
              temp_rho=0.

              IF (dim.EQ.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)

              rho(i) = rho(i) +(dt/2.)*( drho(i)+ temp_rho)

           END IF                               ! End If Inner-Inner-Inner-Loop 2.2.1.1

!          Do Inner-Inner-Inner-Loop 2.2.1.2
           DO d = 1, dim

              v_min(d, i) = vx(d, i)
              vx(d, i) = vx(d, i) + (dt/2.)*dvx(d, i)

           END DO                               ! End Do Inner-Inner-Inner-Loop 2.2.1.2

        END DO                                  ! End Do Inner-Inner-Loop 2.2.1

     END IF                                     ! End If Inner-Loop 2.2

! Definition of variables out of the function vector:

     CALL single_step ( itimestep, dt, ntotal, hsml, mass, x, vx, u, s, rho, p, t, tdsdt,&
                           dx, dvx, du, ds, drho, itype, av)

!    If Inner-Loop 2.3
     IF (itimestep .EQ. 1) THEN

!       Do Inner-Inner-Loop 2.3.1
        DO i=1,ntotal

           temp_u=0.

           IF (dim.EQ.1) temp_u=-nsym*p(i)*vx(1,i)/x(1,i)/rho(i)

           u(i) = u(i) + (dt/2.)*(du(i) + temp_u)

           IF (u(i).LT.0) u(i) = 0.

!          If Inner-Inner-Inner-Loop 2.3.1.1
           IF (.NOT.summation_density ) THEN

              temp_rho=0.

              IF (dim.EQ.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)
              
              rho(i) = rho(i) + (dt/2.)* (drho(i)+temp_rho)

           END IF                               ! End If Inner-Inner-Inner-Loop 2.3.1.1

!          Do Inner-Inner-Inner-Loop 2.3.1.2
           DO d = 1, dim
        
              vx(d, i) = vx(d, i) + (dt/2.) * dvx(d, i) + av(d, i)
              x(d, i) = x(d, i) + dt * vx(d, i)

           END DO                               ! End Do Inner-Inner-Inner-Loop 2.3.1.2

        END DO                                  ! End Do Inner-Inner-Loop 2.3.1

!    Continue If Inner-Loop 2.3
     ELSE

!       Do Inner-Inner-Loop 2.3.2
        DO i=1,ntotal

           temp_u = 0.

           IF (dim.EQ.1) temp_u=-nsym*p(i)*vx(1,i)/x(1,i)/rho(i)

           u(i) = u_min(i) + dt*(du(i)+temp_u)

           IF (u(i).LT.0) u(i) = 0.

!          If Inner-Inner-Inner-Loop 2.3.2.1
           IF (.NOT.summation_density ) THEN

              temp_rho=0.

              IF (dim.EQ.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)

              rho(i) = rho_min(i) + dt*(drho(i)+temp_rho)

           END IF                               ! End If Inner-Inner-Inner-Loop 2.3.2.1

!          Do Inner-Inner-Inner-Loop 2.3.2.2
           DO d = 1, dim

              vx(d, i) = v_min(d, i) + dt * dvx(d, i) + av(d, i)
              x(d, i) = x(d, i) + dt * vx(d, i)

           END DO                               ! End Do Inner-Inner-Inner-Loop 2.3.2.1

        END DO                                  ! End Do Inner-Inner-Loop 2.3.2

     END IF                                     ! End If Inner-Loop 2.3

     time = time + dt

!    If Inner-Loop 2.4
     IF (mod(itimestep,save_step).EQ.0) THEN

        CALL output( x, vx, mass, rho, p, u, c, itype, hsml, ntotal ,dt)          !included dt

     END IF                                     ! End If Inner-Loop 2.4

!    If Inner-Loop 2.5
     IF (mod(itimestep,print_step).EQ.0) THEN

        WRITE(*,*)
        WRITE(*,101)'x','velocity', 'dvx'
        WRITE(*,100)x(1,moni_particle), vx(1,moni_particle), dvx(1,moni_particle)

     END IF                                     ! End If Inner-Loop 2.5

        101 FORMAT(1x,3(2x,a12))
        100 FORMAT(1x,3(2x,e13.6))

  END DO                                        ! End Do Loop 2

nstart=current_ts

END SUBROUTINE time_integration
