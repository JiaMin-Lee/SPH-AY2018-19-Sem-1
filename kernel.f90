SUBROUTINE kernel ( r, dx, hsml, w, dwdx )
! 
! Subroutine to calculate the smoothing kernel wij and its derivatives dwdxij.
!    if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!           = 2, Gauss kernel (Gingold and Monaghan 1981)
!           = 3, Quintic kernel (Morris 1997)
! 
!    r    : Distance between particles i and j                                   [in]
!    dx   : x-, y- and z-distance between i and j                                [in]
!    hsml : Smoothing length                                                     [in]
!    w    : Kernel for all interaction pairs                                    [out]
!    dwdx : Derivative of kernel with respect to x, y and z                     [out]
!
IMPLICIT NONE 

INCLUDE 'param.inc.txt'

! Data dictionary: declare calling parameter types

REAL(KIND=8), INTENT(IN) :: r
REAL(KIND=8), INTENT(IN) :: dx(dim)
REAL(KIND=8), INTENT(IN) :: hsml
REAL(KIND=8), INTENT(OUT) :: w
REAL(KIND=8), INTENT(OUT) :: dwdx(dim)

! Data dictionary: declare local variable types

INTEGER :: i, j, d
REAL(KIND=8) :: q, dw, factor


q = r/hsml
w = 0.E0

! Do Loop 1
  DO d=1,dim

  dwdx(d) = 0.E0

  END DO                                        ! End Do Loop 1

! If Loop 2
  IF (skf.EQ.1) THEN 

!    If Inner-Loop 2.1
     IF (dim.EQ.1) THEN

        factor = 1.E0/hsml

!    Continue If Inner-Loop 2.1
     ELSE IF (dim.EQ.2) THEN

        factor = 15.E0/(7.E0*pi*hsml*hsml)

!    Continue If Inner-Loop 2.1
     ELSE IF (dim.EQ.3) THEN

        factor = 3.E0/(2.E0*pi*hsml*hsml*hsml)

!    Contine If Inner-Loop 2.1
     ELSE
        WRITE (*,*)' >>> Error <<< : Wrong dimension: Dim =',dim
        
        STOP                                      ! Stop program execution of program

     END IF                                     ! End If Loop 2.1

!    If Inner-Loop 2.2
     IF (q.GE.0.AND.q.LE.1.E0) THEN

        w = factor * (2./3. - q*q + q**3 / 2.)

!       Do Inner-Inner-Loop 2.2.1
        DO d = 1, dim

           dwdx(d) = factor * (-2.+3./2.*q)/hsml**2 * dx(d)

        END DO                                  ! End Do Inner-Inner-Loop 2.2.1

!    Continue If Inner-Loop 2.2
     ELSE IF (q.GT.1.E0.AND.q.LE.2) THEN

        w = factor * 1.E0/6.E0 * (2.-q)**3
        
!       Do Inner-Inner-Loop 2.2.2
        DO d = 1, dim

           dwdx(d) =-factor * 1.E0/6.30 * 3.*(2.-q)**2/hsml * (dx(d)/r)

        END DO                                  ! End Do Inner-Inner-Loop 2.2.2

!    Continue If Inner-Loop 2.2
     ELSE

        w=0.

!       Do Inner-Inner-Loop 2.2.3
        DO d= 1, dim

           dwdx(d) = 0.

        END DO                                  ! End Do Inner-Inner-Loop 2.2.3

     END IF                                     ! End If Inner-Loop 2.2

! Continue If Loop 2
  ELSE IF (skf.EQ.2) THEN

     factor = 1.E0 / (hsml**dim * pi**(dim/2.))

!    If Inner-Loop 2.3
     IF (q.GE.0.AND.q.LE.3) THEN

        w = factor * exp(-q*q)

!       Do Inner-Inner-Loop 2.3.1
        DO d = 1, dim

           dwdx(d) = w * ( -2.* dx(d)/hsml/hsml)

        END DO                                  ! End Do Inner-Inner-Loop 2.3.1

!    Continue If Inner-Loop 2.3
     ELSE

        w = 0.

!       Do Inner-Inner-Loop 2.3.2
        DO d = 1, dim

           dwdx(d) = 0.

        END DO                                  ! End Do Inner-Inner-Loop 2.3.2

     END IF                                     ! End If Inner-Loop 2.3

! Continue If Loop 2
  ELSE IF (skf.EQ.3) THEN

!    If Inner-Loop 2.4
     IF (dim.EQ.1) THEN

        factor = 1.E0 / (120.E0*hsml)

!    Continue Inner-Loop 2.4
     ELSE IF (dim.EQ.2) THEN

        factor = 7.e0 / (478.E0*pi*hsml*hsml)

!    Continue Inner-Loop 2.4
     ELSE IF (dim.EQ.3) THEN

        factor = 1.E0 / (120.E0*pi*hsml*hsml*hsml)

!    Continue Inner-Loop 2.4
     ELSE
     
        WRITE (*,*)' >>> Error <<< : Wrong dimension: Dim =', dim
     
        STOP                                    ! Stop program from executing

     END IF                                     ! End If Inner-Loop 2.4

!    If Inner-Loop 2.5
     IF (q.GE.0.AND.q.LE.1) THEN

        w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )

!       Do Inner-Inner-Loop 2.5.1
        DO d = 1, dim

           dwdx(d) = factor * ( (-120 + 120*q - 50*q**2) / hsml**2 * dx(d) )

        END DO                                  ! End Do Inner-Inner-Loop 2.5.1

!    Continue If Inner-Loop 2.5
     ELSE IF (q.GT.1.AND.q.LE.2) THEN

        w = factor * ( (3-q)**5 - 6*(2-q)**5 )

!       Do Inner-Inner-Loop 2.5.2
        DO d = 1, dim

        dwdx(d) = factor * (-5*(3-q)**4 + 30*(2-q)**4) / hsml * (dx(d)/r)

        END DO                                  ! End Do Inner-Inner-Loop 2.5.2

!    Continue If Inner-Loop 2.5
     ELSE IF (q.GT.2.AND.q.LE.3) THEN

        w = factor * (3-q)**5

!       Do Inner-Inner-Loop 2.5.3
        DO d = 1, dim

           dwdx(d) = factor * (-5* (3-q) **4) / hsml * (dx(d)/r)

        END DO                                  ! End Do Inner-Inner-Loop 2.5.3

!    Continue If Inner-Loop 2.5
     ELSE                                      

        w = 0.

!       Do Inner-Inner-Loop 2.5.4
        DO d = 1, dim

           dwdx(d) = 0.

        END DO                                  ! End Do Inner-Inner-Loop 2.5.4

     END IF                                     ! End If Inner-Loop 2.5

  END IF                                        ! End If Loop 2

END SUBROUTINE kernel
