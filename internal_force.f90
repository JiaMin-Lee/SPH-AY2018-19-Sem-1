SUBROUTINE int_force (itimestep,dt,ntotal,hsml,mass,vx,niac,rho,eta,pair_i,pair_j,dwdx,u,itype,x,t,c,p,dvxdt,tdsdt,dedt)
!
! Subroutine to calculate the internal forces on the right hand side
! of the Navier-Stokes equations, i.e. the pressure gradient and the
! gradient of the viscous stress tensor, used by the time integration,
! Moreover the entropy production due to viscous dissipation, tds/dt,
! and the change of internal energy per mass, de/dt, are calculated.
!
! itimestep: Current timestep number                                    [in]
! dt     : Time step                                                    [in]
! ntotal : Number of particles                                          [in]
! hsml   : Smoothing Length                                             [in]
! mass   : Particle masses                                              [in]
! vx     : Velocities of all particles                                  [in]
! niac   : Number of interaction pairs                                  [in]
! rho    : Density                                                      [in]
! eta    : Dynamic viscosity                                            [in]
! pair_i : List of first partner of interaction pair                    [in]
! pair_j : List of second partner of interaction pair                   [in]
! dwdx   : Derivative of kernel with respect to x, y and z              [in]
! itype  : Type of particle (material types)                            [in]
! u      : Particle internal energy                                     [in]
! x      : Particle coordinates                                         [in]
! itype  : Particle type                                                [in]
! t      : Particle temperature                                     [in/out]
! c      : Particle sound speed                                        [out]
! p      : Particle pressure                                           [out]
! dvxdt  : Acceleration with respect to x, y and z                     [out]
! tdsdt  : Production of viscous entropy                               [out]
! dedt   : Change of specific internal energy                          [out]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

INTEGER, INTENT(IN) :: itimestep
INTEGER, INTENT(IN) :: ntotal
INTEGER, INTENT(IN) :: niac
INTEGER, INTENT(IN) :: pair_i(max_interaction)
INTEGER, INTENT(IN) :: pair_j(max_interaction)
INTEGER, INTENT(IN) :: itype(maxn)

REAL(KIND=8), INTENT(IN) :: dt
REAL(KIND=8), INTENT(IN) :: hsml(maxn)
REAL(KIND=8), INTENT(IN) :: mass(maxn)
REAL(KIND=8), INTENT(IN) :: vx(dim,maxn)
REAL(KIND=8), INTENT(IN) :: rho(maxn)
REAL(KIND=8), INTENT(IN) :: eta(maxn)
REAL(KIND=8), INTENT(IN) :: dwdx(dim,max_interaction)
REAL(KIND=8), INTENT(IN) :: u(maxn)
REAL(KIND=8), INTENT(IN) :: x(dim,maxn)
REAL(KIND=8), INTENT(INOUT) :: t(maxn)
REAL(KIND=8), INTENT(OUT) :: c(maxn)
REAL(KIND=8), INTENT(OUT) :: p(maxn)
REAL(KIND=8), INTENT(OUT) :: dvxdt(dim,maxn)
REAL(KIND=8), INTENT(OUT) :: tdsdt(maxn)
REAL(KIND=8), INTENT(OUT) :: dedt(maxn)

! Data dictionary: declare local variable types

INTEGER :: i, j, k, d

REAL(KIND=8) :: dvx(dim), txx(maxn), tyy(maxn), tzz (maxn) , txy(maxn), txz (maxn) , tyz (maxn) , vcc (maxn)
REAL(KIND=8) :: hxx, hyy, hzz, hxy, hxz, hyz, h, hvcc, he, rhoij

! Initialization of shear tensor, velocity divergence, viscous energy, internal energy, acceleration

! Do Loop 1
  DO i=1,ntotal

     txx(i) = 0.e0
     tyy(i) = 0.E0
     tzz(i) = 0.E0
     txy(i) = 0.E0
     txz(i) = 0.E0
     tyz(i) = 0.E0
     vcc(i) = 0.E0
     tdsdt(i) = 0.E0
     dedt(i) = 0.E0

!    Do Inner-Loop 1.1
     DO d=1,dim

        dvxdt(d,i) = 0.E0

     END DO                                     ! End Do Inner-Loop 1.1

  END DO                                        ! End Do Loop 1

! Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc, c

! If Loop 2
  IF (visc) THEN

!    Do Inner-Loop 2.1
     DO k=1,niac

        i = pair_i(k)
        j = pair_j(k)

!       Do Inner-Inner-Loop 2.1.1
        DO d=1,dim

           dvx(d) = vx(d,j) - vx(d,i)

        END DO                                  ! End Do Inner-Inner-Loop 2.1.1

!       If Inner-Inner-Loop 2.1.2 (Block 1)
        IF (dim.EQ.1) THEN

           hxx = 2.E0*dvx(1)*dwdx(1,k)

!       Continue If Inner-Inner-Loop 2.1.2 (Block 2)
        ELSE IF (dim.EQ.2) THEN

           hxx = 2.E0*dvx(1)*dwdx(1,k) - dvx(2)*dwdx(2,k)
           hxy = dvx(1)*dwdx(2,k) + dvx(2) *dwdx(1,k)
           hyy = 2.E0*dvx(2)*dwdx(2,k) - dvx(1)*dwdx(1,k)

!       Continue If Inner-Inner-Loop 2.1.2 (Block 3)
        ELSE IF (dim.EQ.3) THEN

           hxx = 2.E0*dvx(1)*dwdx(1,k) - dvx(2)*dwdx(2,k) - dvx(3)*dwdx(3,k)
           hxy = dvx(1)*dwdx(2,k) + dvx(2)*dwdx(1,k)
           hxz = dvx(1)*dwdx(3,k) + dvx(3)*dwdx(1,k)
           hyy = 2.E0*dvx(2)*dwdx(2,k) - dvx(1)*dwdx(1,k) - dvx(3)*dwdx(3,k)
           hyz = dvx(2)*dwdx(3,k) + dvx(3)*dwdx(2,k)
           hzz = 2.E0*dvx(3)*dwdx(3,k) - dvx(1)*dwdx(1,k) - dvx(2)*dwdx(2,k)

        END IF                                   ! End If Inner-Inner-Loop 2.1.2

        hxx = 2.E0/3.E0*hxx
        hyy = 2.E0/3.E0*hyy
        hzz = 2.E0/3.E0*hzz

!       If Inner-Inner-Loop 2.1.3 (Block 1)
        IF (dim.EQ.1) THEN

           txx(i) = txx(i) + mass(j)*hxx/rho(j)
           txx(j) = txx(j) + mass(i)*hxx/rho(i)
   
!       Continue If Inner-Inner-Loop 2.1.3 (Block 2)
        ELSE IF (dim.EQ.2) THEN

           txx(i) = txx(i) + mass(j)*hxx/rho(j)
           txx(j) = txx(j) + mass(i)*hxx/rho(i)
           txy(i) = txy(i) + mass(j)*hxy/rho(j)
           txy(j) = txy(j) + mass(i)*hxy/rho(i)
           tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
           tyy(j) = tyy(j) + mass(i)*hyy/rho(i)

!       Continue If Inner-Inner-Loop 2.1.3 (Block 3)
        ELSE IF (dim.EQ.3) THEN

           txx(i) = txx(i) + mass(j)*hxx/rho(j)
           txx(j) = txx(j) + mass(i)*hxx/rho(i)
           txy(i) = txy(i) + mass(j)*hxy/rho(j)
           txy(j) = txy(j) + mass(i)*hxy/rho(i)
           txz(i) = txz(i) + mass(j)*hxz/rho(j)
           txz(j) = txz(j) + mass(i)*hxz/rho(i)
           tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
           tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
           tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
           tyz(j) = tyz(j) + mass(i)*hyz/rho(i)
           tzz(i) = tzz(i) + mass(j)*hzz/rho(j)
           tzz(j) = tzz(j) + mass(i)*hzz/rho(i)
 
        END IF                                     ! End If Inner-Inner-Loop 2.1.3

! Calculate SPH sum for vc,c = dvx/dx + dvy/dy + dvz/dz:

        hvcc = 0.

!       Do  Inner-Inner-Loop 2.1.4
        DO d=1,dim

           hvcc = hvcc + dvx(d)*dwdx(d,k)

        END DO                                     ! End Do Inner-Inner-Loop 2.1.4

        vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
        vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)

     END DO                                        ! End Do Inner-Loop 2.1

  END IF                                           ! End If Loop 2

! Do Loop 3
  DO i=1,ntotal

!    Viscous entropy Tds/dt - 1/2 eta/rho Tab Tab

!    If Inner-Loop 3.1
     IF (visc) THEN

!       If Inner-Inner-Loop 3.1.1
        IF (dim.EQ.1) THEN

           tdsdt(i) = txx(i)*txx(i)

!       Continue If Inner-Inner-Loop 3.1.1
        ELSE IF (dim.EQ.2) THEN

           tdsdt(i) = txx(i)*txx(i) + 2.E0*txy(i)*txy(i) + tyy(i)*tyy(i)

!       Continue If Inner-Inner-Loop 3.1.1
        ELSE IF (dim.EQ.3) THEN

           tdsdt(i) = txx(i)*txx(i) + 2.E0*txy(i)*txy(i) + 2.E0*txz(i)*txz(i) + tyy(i)*tyy(i) + 2.E0*tyz(i)*tyz(i)  + tzz(i)*tzz(i)

        END IF                                      ! End If Inner-Inner-Loop 3.1.1

        tdsdt(i) = 0.5e0*eta(i)/rho(i)*tdsdt(i)

     END IF                                         ! End If Inner-Loop 3.1

! Pressure from equation of state

!    If Inner-Loop 3.2
     IF (abs(itype(i)).EQ.1) THEN

        call p_gas(rho(i), u(i), p(i),c(i))

     ELSE IF (abs(itype(i)).EQ.2) THEN

        call p_art_water(rho(i), p(i), c(i))

     END IF                                         ! End If Inner-Loop 3.2

  END DO                                            ! End Do Loop 3

! Calculate SPH sum for pressure force -p,a/rho
! and viscous force (eta Tab),b/rho
! and the internal energy change de/dt due to -p/rho vc,c

! Do Loop 4
  DO k=1,niac

     i = pair_i(k)
     j = pair_j(k)
     he = 0.E0

!    For SPH algorithm 1

     rhoij = 1.E0/(rho(i)*rho(j))

!    If Inner-Loop 4.1
     IF (pa_sph.EQ.1) THEN

!       Do Inner-Inner-Loop 4.1.1
        DO d=1,dim

!          Pressure part

           h = -(p(i) + p(j))*dwdx(d,k)
           he = he + (vx(d,j) - vx(d,i))*h

!          Viscous force

!          If Inner-Inner-Inner-Loop 4.1.1.1
           IF (visc) THEN

!             If Inner-Inner-Inner-Inner-Loop 4.1.1.1.1
              IF (d.EQ.1) THEN

!                x-coordinate of acceleration

                 h = h + (eta(i)*txx(i) + eta(j)*txx(j))*dwdx(1,k)

!                If Inner-Inner-Inner-Inner-Inner-Loop 4.1.1.1.1.1
                 IF (dim.GE.2) THEN

                    h = h + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(2,k)

!                   If Inner-Inner-Inner-Inner-Inner-Inner-Loop 4.1.1.1.1.1.1
                    IF (dim.EQ.3) THEN

                        h = h + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(3,k)

                    END IF                                      ! End If 4.1.1.1.1.1.1

                 END IF                                         ! End If 4.1.1.1.1.1

!             Continue If Inner-Inner-Inner-Inner-Loop 4.1.1.1.1
              ELSE IF (d.EQ.2) THEN

!                y-coordinate of acceleration

                 h = h + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(1,k) + (eta(i)*tyy(i) + eta(j)*tyy(j))*dwdx(2,k)

!                If Inner-Inner-Inner-Inner-Inner-Loop 4.1.1.1.1.2
                 IF (dim.EQ.3) THEN

                    h = h + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(3,k)

                 END IF                                         ! End If 4.1.1.1.1.2

!             Continue If Inner-Inner-Inner-Inner-Loop 4.1.1.1.1
              ELSE IF (d.EQ.3) THEN

!                z-coordinate of acceleration

                 h = h + (eta(i)*txz(i)+eta(j)*txz(j))*dwdx(1,k)&
                       + (eta(i)*tyz(i)+eta(j)*tyz(j))*dwdx(2,k)&
                       + (eta(i)*tzz(i)+eta(j)*tzz(j))*dwdx(3,k)

              END IF                                            ! End If 4.1.1.1.1
        
           END IF                                               ! End If 4.1.1.1

           h = h*rhoij
           dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
           dvxdt(d,j) = dvxdt(d,j) - mass(i)*h

        END DO                                                  ! End Do 4.1.1

        he = he*rhoij
        dedt(i) = dedt(i) + mass(j)*he
        dedt(j) = dedt(j) + mass(i)*he

!    For SPH algorithm 2

!    Continue If Inner-Loop 4.1
     ELSE IF (pa_sph.EQ.2) THEN

!       Do Inner-Inner-Loop 4.1.2
        DO d=1,dim

           h = - (p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(d,k)
           he = he + (vx(d,j) - vx(d,i))*h

!          Viscous force

!          If Inner-Inner-Inner-Loop 4.1.2.1
           IF (visc) THEN

!             If Inner-Inner-Inner-Inner-Loop 4.1.2.1.1
              IF (d.EQ.1) THEN

!                 x-coordinate of acceleration

                  h = h + (eta(i)*txx(i)/rho(i)**2 + eta(j)*txx(j)/rho(j)**2)*dwdx(1,k)

!                 If Inner-Inner-Inner-Inner-Inner-Loop 4.1.2.1.1.1
                  IF (dim.GE.2) THEN

                     h = h + (eta(i)*txy(i)/rho(i)**2 + eta(j)*txy(j)/rho(j)**2)*dwdx(2,k)

!                    If Inner-Inner-Inner-Inner-Inner-Inner-Loop 4.1.2.1.1.1.1
                     IF (dim.EQ.3) THEN

                        h = h + (eta(i)*txz(i)/rho(i)**2 + eta(j)*txz(j)/rho(j)**2)*dwdx(3,k)
   
                     END IF                                                ! End If 4.1.2.1.1.1.1

                  END IF                                                   ! End If 4.1.2.1.1.1 

!             Continue If 4.1.2.1.1
              ELSE IF (d.EQ.2) THEN

!                 y-coordinate of acceleration

                  h= h + (eta(i)*txy(i)/rho(i)**2 + eta(j)*txy(j)/rho(j)**2)*dwdx(1,k)&
                       + (eta(i)*tyy(i)/rho(i)**2 + eta(j)*tyy(j)/rho(j)**2)*dwdx(2,k)

!                 If Inner-Inner-Inner-Inner-Inner-Loop 4.1.2.1.1.2
                  IF (dim.EQ.3) THEN

                     h = h + (eta(i)*tyz(i)/rho(i)**2 + eta(j)*tyz(j)/rho(j)**2)*dwdx(3,k)

                  END IF                                                  ! End If 4.1.2.1.1.2

!             Continue If 4.1.2.1.1
              ELSE IF (d.EQ.3) THEN

!                 z-coordinate of acceleration

                  h = h + (eta(i)*txz(i)/rho(i)**2 + eta(j)*txz(j)/rho(j)**2)*dwdx(1,k)&
                        + (eta(i)*tyz(i)/rho(i)**2 + eta(j)*tyz(j)/rho(j)**2)*dwdx(2,k)&
                        + (eta(i)*tzz(i)/rho(i)**2 + eta(j)*tzz(j)/rho(j)**2)*dwdx(3,k)

              END IF                                                      ! End If 4.1.2.1.1

           END IF                                                         ! End If 4.1.2.1

           dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
           dvxdt(d,j) = dvxdt(d,j) - mass(i)*h

        END DO                                                            ! End Do 4.1.2

        dedt(i) = dedt(i) + mass(j)*he
        dedt(j) = dedt(j) + mass(i)*he

     END IF                                                               ! End If 4.1

  END DO                                                                  ! End Do Loop 4

! Change of specific internal energy de/dt = T ds/dt - p/rho vc,c:

! Do Loop 5
  DO i=1,ntotal

     dedt(i) = tdsdt(i) + 0.5E0*dedt(i)

  END DO                                                        ! End Do Loop 5

END SUBROUTINE int_force
