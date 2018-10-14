SUBROUTINE direct_find ( itimestep, ntotal , hsml, x, niac, pair_i, pair_j, w, dwdx, countiac )
!
! Subroutine to calculate the smoothing funciton for each particle and the interaction parameters used by the 
! SPH algorithm. Interaction pairs are determined by directly comparing the particle distance with the
! corresponding smoothing length.

! itimestep : Current time step                                         [in]
! ntotal    : Number of particles                                       [in]
! hsml      : Smoothing Length                                          [in]
! x         : Coordinates of all particles                              [in]
! niac      : Number of interaction pairs                              [out]
! pair_i    : List of first partner of interaction pair                [out]
! pair_j    : List of second partner of interaction pair               [out]
! w         : Kernel for all interaction pairs                         [out]
! dwdx      : Derivative of kernel with respect to x, y and z          [out]
! countiac  : Number of neighboring particles                          [out]

IMPLICIT NONE

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types
INTEGER, INTENT(IN) :: itimestep
INTEGER, INTENT(IN) :: ntotal
INTEGER, INTENT(OUT) :: niac
INTEGER, INTENT(OUT) :: pair_i(max_interaction)
INTEGER, INTENT(OUT) :: pair_j(max_interaction)
INTEGER, INTENT(OUT) :: countiac(maxn)

REAL(KIND=8), INTENT(IN) :: hsml(maxn)
REAL(KIND=8), INTENT(IN) :: x(dim,maxn)
REAL(KIND=8), INTENT(OUT) :: w(max_interaction)
REAL(KIND=8), INTENT(OUT) :: dwdx(dim,max_interaction)

! Data dictionary: declare local variable types
  
INTEGER :: i, j, d, sumiac, maxiac, miniac, noiac, maxp, minp, scale_k

REAL(KIND=8) :: dxiac(dim), driac, r, mhsml, tdwdx(dim)

! If Loop 1
  IF (skf.EQ.1) THEN
     scale_k = 2
  ELSE IF (skf.EQ.2) THEN
     scale_k = 3
  ELSE IF (skf.EQ.3) THEN
     scale_k = 3
  END IF                                        ! End If Loop 1

! Do Loop 2
  DO i=1,ntotal
     countiac(i) = 0
  END DO                                        ! End Do Loop 2

niac = 0

! Do Loop 3
  DO i=1,ntotal-1

!    Do Inner-Loop 3.1
     DO j = i+1, ntotal
        dxiac(1) = x(1,i) - x(1,j)
        driac = dxiac(1)*dxiac(1)

!       Do Inner-inner-Loop 3.1.1    
        DO d=2,dim
           dxiac(d) = x(d,i) - x(d,j)
           driac = driac + dxiac(d)*dxiac(d)
        END DO                                  ! End Do Inner-Inner-Loop 3.1.1

     mhsml = (hsml(i)+hsml(j))/2.

!    If Inner-Loop 3.2
     IF (sqrt(driac).LT.scale_k*mhsml) THEN

!       If Inner-Inner Loop 3.2.2
        IF (niac.LT.max_interaction) THEN

!          Neighboring pair list, and totalinteraction number and the interaction number for each particle
           niac = niac + 1
           pair_i(niac) = i
           pair_j(niac) = j
           r = sqrt(driac)
           countiac(i) = countiac(i) + 1
           countiac(j) = countiac(j) + 1

!          Kernel and derivations of kernel
!          Call subroutine kernel
           CALL kernel(r, dxiac, mhsml, w(niac), tdwdx)

!          Do Inner-Inner-Inner-Loop 3.2.2.1
           DO d=1,dim
              dwdx(d,niac) = tdwdx(d)
           END DO                               ! End Do Inner-Inner-Inner-Loop 3.2.2.1

        ELSE
           WRITE (*,*) ' >>> ERROR <<< : Too many interactions'        
           STOP                                 ! Stop program
        END IF                                  ! End If Inner-Inner-Loop 3.2.2

     END IF                                     ! End If Inner-Loop 3.2
   
     END DO                                     ! End Do Inner-Loop 3.1

  END DO                                        ! End Do Loop 3

! Statistics for the interaction

  sumiac = 0
  maxiac = 0
  miniac = 1000
  noiac = 0

! Do Loop 4
  DO i=1,ntotal
     sumiac = sumiac + countiac(i)

!    If Inner-Loop 4.1
     IF (countiac(i).GT.maxiac) THEN
        maxiac = countiac(i)
        maxp = i
     END IF                                     ! End If Inner-Loop 4.1

!    If Inner-Loop 4.2
     IF (countiac(i).LT.miniac) THEN
        miniac = countiac(i)
        minp = i
     END IF                                     ! End If Inner-Loop 4.2

!    If Inner-Loop 4.3                         
     IF (countiac(i).EQ.0) THEN
        noiac = noiac + 1
     END IF                                     ! End If Inner-Loop 4.3
  
  END DO                                        ! End Do Loop 4

! If Loop 5
  IF (mod(itimestep,print_step).EQ.0) THEN

!    If Loop 5.1
     IF (int_stat) THEN
        WRITE (*,*) ' >> Statistics: interactions per particle:'
        WRITE (*,*) '**** Particle:',maxp, ' maximal interactions:',maxiac
        WRITE (*,*) '**** Particle:',minp, ' minimal interactions:',miniac
        WRITE (*,*) '**** Average :',real(sumiac)/real(ntotal)
        WRITE (*,*) '**** Total pairs : ',niac
        WRITE (*,*) '**** Particles with no interactions:',noiac
     END IF                                     ! End If Loop 5.1

  END IF                                        ! End If Loop 5

END SUBROUTINE direct_find
