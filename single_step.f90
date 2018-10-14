SUBROUTINE single_step ( itimestep, dt, ntotal, hsml, mass, x, vx, u, s, rho, p, t, tdsdt, dx, dvx, du, ds, drho,&
                         itype, av )
!
! Subroutine to determine the right hand side of a differential
! equation in a single step for performing time integration
! In this routine and its subroutines the SPH algorithms are performed.
!
! itimestep : Current timestep number                                     [in]
! dt        : Timestep                                                    [in]
! ntotal    : Number of particles                                         [in]
! hsml      : Smoothing Length                                            [in]
! mass      : Particle masses                                             [in]
! x         : Particle position                                           [in]
! vx        : Particle velocity                                           [in]
! u         : Particle internal energy                                    [in]
! s         : Particle entropy (not used here)                            [in]
! rho       : Density                                                 [in/out]
! p         : Pressure                                                   [out]
! t         : Temperature                                             [in/out]
! tdsdt     : Production of viscous entropy t*ds/dt                      [out]
! dx        : dx = vx = dx/dt                                            [out]
! dvx       : dvx = dvx/dt, force per unit mass                          [out]
! du        : du = du/dt                                                 [out]
! ds        : ds = ds/dt                                                 [out]
! drho      : drho = drh,o/dt                                            [out]
! itype     : Type of particle                                            [in]
! av        : Monaghan average velocity                                  [out]
!
IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

INTEGER, INTENT(IN) :: itimestep
INTEGER, INTENT(IN) :: ntotal
INTEGER, INTENT(IN) :: itype(maxn)

REAL(KIND=8), INTENT(IN) :: dt
REAL(KIND=8), INTENT(IN) :: hsml(maxn)
REAL(KIND=8), INTENT(IN) :: mass(maxn)
REAL(KIND=8), INTENT(IN) :: x(dim,maxn)
REAL(KIND=8), INTENT(IN) :: vx(dim,maxn)
REAL(KIND=8), INTENT(IN) :: u(maxn)
REAL(KIND=8), INTENT(IN) :: s(maxn)
REAL(KIND=8), INTENT(INOUT) :: rho(maxn)
REAL(KIND=8), INTENT(OUT) :: p(maxn)
REAL(KIND=8), INTENT(INOUT) :: t(maxn)
REAL(KIND=8), INTENT(OUT) :: tdsdt(maxn)
REAL(KIND=8), INTENT(OUT) :: dx(dim,maxn)
REAL(KIND=8), INTENT(OUT) :: dvx(dim,maxn)
REAL(KIND=8), INTENT(OUT) :: du(maxn)
REAL(KIND=8), INTENT(OUT) :: ds(maxn)
REAL(KIND=8), INTENT(OUT) :: drho(maxn)
REAL(KIND=8), INTENT(OUT) :: av(dim, maxn)

! Data dictionary: declare local variable types

INTEGER :: i, d, nvirt, niac, pair_i(max_interaction) , pair_j(max_interaction), ns(maxn)

REAL(KIND=8) :: w(max_interaction), dwdx(dim,max_interaction), indvxdt(dim,maxn), exdvxdt(dim,maxn),ardvxdt(dim,maxn)
REAL(KIND=8) :: avdudt(maxn), ahdudt(maxn), c(maxn), eta(maxn)

! Do Loop 1
  DO i=1,ntotal

     avdudt(i) = 0.
     ahdudt(i) = 0.

!    Do Inner-Loop 1.1
     DO d=1,dim
        indvxdt(d,i) = 0.
        ardvxdt(d,i) = 0.
        exdvxdt(d,i) = 0.

     END DO                                     ! End Do Inner-Loop 1.1

  END DO                                        ! End Do Loop 1

! Positions of virtual (boundary) particles:

  nvirt = 0

! If Loop 2
  IF (virtual_part) THEN

     CALL virt_part ( itimestep, ntotal, nvirt, hsml, mass, x, vx, rho, u, p, itype )

  END IF                                        ! End If Loop 2

! Interaction parameters, calculating neighboring particles and optimzing smoothing length

! If Loop 3
  IF (nnps.EQ.1) THEN

     CALL direct_find ( itimestep, ntotal+nvirt, hsml, x, niac, pair_i, pair_j, w, dwdx, ns )

! Continue If Loop 3
  ELSE IF (nnps.EQ.2) THEN

     CALL link_list ( itimestep, ntotal+nvirt, hsml(1), x, niac, pair_i, pair_j, w, dwdx, ns )

! Continue If Loop 3
  ELSE IF (nnps.EQ.3) THEN

     CALL tree_search ( itimestep, ntotal+nvirt, hsml, x, niac,pair_i, pair_j, w, dwdx, ns )

  END IF                                        ! End If Loop 3

! Density approximation or change rate

! If Loop 4
  IF (summation_density) THEN

     CALL sum_density ( ntotal+nvirt, hsml, mass, niac, pair_i, pair_j, w, itype, rho )

! Continue Loop 4
  ELSE

     CALL con_density ( ntotal+nvirt, mass, niac, pair_i, pair_j, dwdx, vx, itype, x, rho, drho )
  
  END IF                                        ! End If Loop 4

! Dynamic viscosity:

  IF (visc) CALL  viscosity ( ntotal+nvirt, itype, x, rho, eta )

! Internal forces:

  CALL int_force ( itimestep, dt, ntotal+nvirt, hsml, mass, vx, niac, rho, eta, pair_i, pair_j, dwdx, u, itype,&
                  x, t, c, p, indvxdt, tdsdt, du )

! Artificial viscosity:

  IF (visc_artificial) CALL art_visc ( ntotal+nvirt, hsml, mass, x, vx, niac, rho, c, pair_i, pair_j, w, dwdx,&
                                       ardvxdt, avdudt )

! External forces:

  IF (ex_force) CALL ext_force ( ntotal+nvirt, mass, x, niac, pair_i, pair_j, itype, hsml, exdvxdt)

! Calculating Che neighboring particles and undating HSML

  IF (sle.NE.0) CALL h_upgrade ( dt, ntotal, mass, vx, rho, niac, pair_i, pair_j, dwdx, hsml )

  IF (heat_artificial) CALL art_heat ( ntotal+nvirt, hsml, mass, x, vx, niac, rho, u, c, pair_i, pair_j, w,&
                                       dwdx, ahdudt)

! Calculating average velocity of each partile for avoiding penetration

  IF (average_velocity) CALL av_vel ( ntotal, mass, niac, pair_i, pair_j, w, vx, rho, av)

! Convert velocity, force, and energy to f and dfdt

! Do Loop 5
  DO i=1,ntotal

!    Do Inner-Loop 5.1
     DO d=1,dim

        dvx(d,i) = indvxdt(d,i) + exdvxdt(d,i) + ardvxdt(d,i)

     END DO                                     ! End Do Inner-Loop 5.1

     du(i) = du(i) + avdudt(i) + ahdudt(i)

  END DO                                        ! End Do Loop 5

! If Loop 6
  IF (mod ( itimestep, print_step ).EQ.0) THEN

     WRITE(*,*)
     WRITE(*,*) '**** Information for particle ****', moni_particle
     WRITE(*,101)'internal a ','artifical a=', 'external a ' , 'total a '
     WRITE(*,100)indvxdt(1,moni_particle),ardvxdt(1,moni_particle),exdvxdt(1,moni_particle),dvx(1,moni_particle)

  END IF                                        ! End If Loop 6

101 FORMAT(1x,4(2x,A12))
100 FORMAT(1x,4(2x,E13.6))

END SUBROUTINE single_step
