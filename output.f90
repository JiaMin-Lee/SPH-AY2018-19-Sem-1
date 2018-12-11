SUBROUTINE output( x, vx, mass, rho, p, u, c, itype, hsml, ntotal, maxtimestep, numx, numy, dt )
!
! Subroutine for saving particle information to external disk file
!
! x      : coordinates of particles                             [in]
! vx     : velocities of particles                              [in]
! mass   : mass of particles                                    [in]
! rho    : densities of particles                               [in]
! p      : pressure of particles                                [in]
! u      : internal energy of particles                         [in]
! c      : sound velocity of particles                          [in]
! itype  : types of particles                                   [in]
! hsml   : smoothing lengths of particles                       [in]
! ntotal : total real  particle number                          [in]
!

use loguns

IMPLICIT NONE

INCLUDE 'param.inc.txt'


! Data dictionary: declare calling parameter types

INTEGER, INTENT(IN) :: itype(maxn)
INTEGER, INTENT(IN) :: ntotal
!INTEGER, INTENT(IN) :: nvirt                    ! INCLUDED nvirt
INTEGER, INTENT(IN) :: maxtimestep               ! INCLUDED maxtimestep
INTEGER, INTENT(IN) :: numx, numy                ! INCLUDED num

REAL(KIND=8), INTENT(IN) :: x(dim, maxn)
REAL(KIND=8), INTENT(IN) :: vx(dim, maxn)
REAL(KIND=8), INTENT(IN) :: mass(maxn)
REAL(KIND=8), INTENT(IN) :: rho(maxn)
REAL(KIND=8), INTENT(IN) :: p(maxn)
REAL(KIND=8), INTENT(IN) :: u(maxn)
REAL(KIND=8), INTENT(IN) :: c(maxn)
REAL(KIND=8), INTENT(IN) :: hsml(maxn)
REAL(KIND=8), INTENT(IN) :: dt                           ! INCLUDED
!REAL(KIND=8), INTENT(IN) :: gamma                       ! INCLUDED

CHARACTER(len=30) :: fxv, fstate, fother, fmeshxv                         ! INCLUDED

!!!!! From Daniel Price
! data dictionary
character(len=len(rootname)+10) :: dumpfile                     ! INCLUDED
real :: gamma, hfact                                            ! INCLUDED
integer :: nprint,ncolumns, ndimV                               ! INCLUDED
integer :: ierr,iformat                                         ! INCLUDED
integer :: imhd, idust, igravity, idumpghost
character(len=12) :: geom


!!!!!!!!!!! From Lius
! Data dictionary: declare local variable types

INTEGER :: i, d, npart                       

!!!!!!! FROM JIAMIN
! Data dictionary: interpolating to mesh 
INTEGER :: mi, miy, mix, n              
REAL(KIND=8) :: x1, y1, lx, ly                         
REAL(KIND=8) :: sum_vxx(maxm), sum_vxy(maxm), avg_vxx(maxm), avg_vxy(maxm)  
REAL(KIND=8) :: mistart(dim,maxm), miend(dim,maxm), mimid(dim,maxm)              

x1 = 1.E-3
y1 = 1.E-3
lx = x1/numx
ly = y1/numy

 write(fxv,"('../data/f_xv',i6.6,'.dat')") maxtimestep
 write(fstate,"('../data/f_state',i6.6,'.dat')") maxtimestep
 write(fother,"('../data/f_other',i6.6,'.dat')") maxtimestep
!write(dumpfile,"('sph_',i5.5,'.dat')") ifile


OPEN(1, file=fxv, status='replace')
OPEN(2, file=fstate, status='replace')
OPEN(3, file=fother, status='replace')


!OPEN(1,file="../data/f_xv.dat") 
!OPEN(2,file="../data/f_state.dat")
!OPEN(3,file="../data/f_other.dat")


!WRITE(1,*) ntotal

! Do Loop 1
  DO i = 1, ntotal

     WRITE(1,1001) i, (x(d, i), d=1,dim), (vx(d, i), d = 1, dim)
     WRITE(2,1002) i, mass(i), rho(i), p(i), u(i)
     WRITE(3,1003) i, itype(i), hsml(i)
  
  END DO                                ! End Do Loop 1

1001 format(1x, I6, 6(2x, E15.8))
1002 format(1x, I6, 7(2x, E15.8))
1003 format(1x, I6, 2x, I4, 2x, E15.8)
! 1001 format(1x, I6, 6(2x, E14.8))
! 1002 format(1x, I6, 7(2x, E14.8))
! 1003 format(1x, I6, 2x, I4, 2x, E14.8)

!CLOSE(1)
!CLOSE(2)
!CLOSE(3)

!!!!!!!! Write data to num x num mesh
! If Loop A1
  IF ((numx.GT.0).AND.(numy.GT.0)) THEN

     WRITE(fmeshxv,"('../data/fmesh_xv',i6.6,'.dat')") maxtimestep

     OPEN(4, file=fmeshxv, status='replace')
   
!    Do Loop A1.1
     DO mix = 1, numx
        
        DO miy = 1, numy   

           mi = miy+(mix-1)*numx

           mistart(1, mi) = (mix-1)*lx
           miend(1,mi)= mix*lx
           mistart(2,mi) = (miy-1)*ly
           miend(2,mi) = miy*ly

           mimid(1,mi)=(mistart(1,mi)+miend(1,mi))/2
           mimid(2,mi)=(mistart(2,mi)+miend(2,mi))/2

           n = 0
           sum_vxx(mi)=0
           sum_vxy(mi)=0

            !READ(1,*)  i, (x(d, i), d=1,dim), (vx(d, i), d = 1, dim)

           DO i=1,ntotal

              IF (x(1,i).GE.mistart(1,mi).AND.x(1,i).LE.miend(1,mi)) THEN
                
                 IF (x(2,i).GE.mistart(2,mi).AND.x(2,i).LE.miend(2,mi)) THEN

                     n=n+1
                     sum_vxx(mi) = sum_vxx(mi)+vx(1,i)
                     sum_vxy(mi) = sum_vxy(mi)+vx(2,i)
                 END IF

              END IF

           END DO     
        
           IF (n.GT.0) THEN 
              avg_vxx(mi) = sum_vxx(mi)/n
              avg_vxy(mi) = sum_vxy(mi)/n
          
           ELSE 
              avg_vxx(mi) = 0 
              avg_vxy(mi) = 0

           END IF

                    
           !WRITE(4,1004) mi, n, (mimid(d, mi), d = 1, dim), avg_vxx(mi), avg_vxy(mi), mistart(1, mi), miend(1,mi), mistart(2,mi), miend(2,mi)
           WRITE(4,1004) mi, n, (mimid(d, mi), d = 1, dim), avg_vxx(mi), avg_vxy(mi)
           1004 FORMAT(1x, I6, 2x, I6,  4(2x, E15.8))
           !1004 FORMAT(1x, I6, 2x, I6, 4(2x, E15.8), 2(E11.4' , 'E11.4))

          

        END DO   

     END DO
    
  ELSE 

         WRITE(*,*) '**********************************************************'
         WRITE(*,*) 'Interpolating to mesh is not selected.'
  
  END IF

!!!!!!!!!!!!! From Daniel Price

imhd = 0 
idust = 0
igravity = 0
idumpghost = 1
geom = 'cartesian'
ndimV = 3
gamma = 1.4
hfact = 2.0       !fac from liu's

! from outputND_mhd.f90
!--create new dumpfile if rootname contains numbers
!

 ifile = ifile + 1
 write(dumpfile,"('../splashdata/sph_',i5.5,'.dat')") ifile

npart = ntotal

! from readwrite_dumps.f90

if (idumpghost.eq.1) then
    nprint = ntotal
 else
    nprint = npart
    nprint = ntotal          
 endif
!
!--open dumpfile
!
 open(unit=idatfile,file=dumpfile,status='replace',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'error: can''t create new dumpfile ',trim(dumpfile)
    stop
 endif

!
!--write timestep header to data file
!
 ncolumns = dim + 2*ndimV + 4
 if (imhd.ne.0) then
    ncolumns = ncolumns + 8 + 2*ndimV ! number of columns
    iformat = 2
    if (imhd.lt.0) ncolumns = ncolumns + ndimV
 else
    ncolumns = ncolumns + 5
    iformat = 1
    if (igravity.ne.0) ncolumns = ncolumns + 1
!    if (allocated(del2v)) ncolumns = ncolumns + 1
 endif
 if (geom(1:4).ne.'cart') then
    ncolumns = ncolumns + 2 + ndimV
    iformat = iformat + 2
 endif
 if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
    iformat = 5
    ncolumns = ncolumns + ndimV + 1
 endif

 write(idatfile,iostat=ierr) dt,npart,nprint,gamma,hfact,dim,ndimV, &
      ncolumns,iformat,x_mingeom, y_mingeom, z_mingeom, x_maxgeom, y_maxgeom, &
      z_maxgeom, len(geom),geom
!removed ibound
!changed xmin(1:dim), xmax(1:dim)
 if (ierr /= 0) then
    write(iprint,*) '*** error writing timestep header to dumpfile',trim(dumpfile)
 endif
!
!--write the data (primitive variables) to the data file
!  data is written in two blocks :

! MHD variables are written after the hydro ones in 
! each case
!
  !--essential variables
  do i=1,dim
     write(idatfile) x(i,1:nprint)
  enddo
  do i=1,dim
     write(idatfile) vx(i,1:nprint)
  enddo
  write(idatfile) hsml(1:nprint)
!  if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
     write(idatfile) rho(1:nprint)
!  else
!     write(idatfile) dens(1:nprint)
!  endif
  write(idatfile) u(1:nprint)
  write(idatfile) mass(1:nprint)

!  if (imhd.ne.0) then
!     do i=1,3
!        write(idatfile) alpha(i,1:nprint)
!     enddo
!     do i=1,ndimV
!        write(idatfile) Bfield(i,1:nprint)
!     enddo
!     write(idatfile) psi(1:nprint)
!     !--info only
!     write(idatfile) pr(1:nprint)
!     write(idatfile) -drhodt(1:nprint)/rho(1:nprint)
!     write(idatfile) divB(1:nprint)
!     do i=1,ndimV
!        write(idatfile) curlB(i,1:nprint)
!     enddo
!     write(idatfile) gradh(1:nprint)
!     do i=1,ndimV
!        write(idatfile) force(i,1:nprint)
!     enddo
!     if (imhd.lt.0) then
!        do i=1,ndimV
!           write(idatfile) Bevol(i,1:nprint)
!        enddo
!        write(idatfile) Bconst(:)
!     endif
!  else
!     do i=1,2
!        write(idatfile) alpha(i,1:nprint)
!     enddo
     !--info only
     write(idatfile) p(1:nprint)
!     write(idatfile) -drhodt(1:nprint)/rho(1:nprint)
!     write(idatfile) dhsml(1:nprint)
!     do i=1,ndimV
!        write(idatfile) force(i,1:nprint)
!     enddo
!     if (igravity.ne.0) write(idatfile) poten(1:nprint)
!     if (allocated(del2v)) write(idatfile) del2v(1:nprint)
!  endif
!  if (geom(1:4).ne.'cart') then
!     write(idatfile) rho(1:nprint)
!     write(idatfile) sqrtg(1:nprint)
!     do i=1,ndimV
!        write(idatfile) pmom(i,1:nprint)
!     enddo
!  endif
!  if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
!     write(idatfile) dustfrac(1:nprint)
!     do i=1,ndimV
!        write(idatfile) deltav(i,1:nprint)
!    enddo
!  endif
  write(idatfile) itype(1:nprint)

 close(unit=idatfile)

! return
!end subroutine write_dump

CLOSE(1)
CLOSE(2)
CLOSE(3)
CLOSE(4)

END SUBROUTINE output
