PROGRAM Integration

!***********************************************************************
!----  This program calculates a defined integralusing differents   ----
!----  methods: Trapezoid, Simpson, Romberg and Gauss-Legendreren   ----
!----                                                               ----
!----                         ~VCastor 2021                         ----
!***********************************************************************
! There is a pdf with documentation about the math theory necessary to
! understand properly the code. explicit_explanations.pdf
!***********************************************************************

!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER                             :: d, n, i
REAL*8                              :: lowl, upl, e, er, s
REAL*8                              :: aux1, aux2, integra
REAL*8, PARAMETER                   :: threshold = 1.d-8
REAL*8, DIMENSION(:), ALLOCATABLE   :: t, w
REAL*8, DIMENSION(:,:), ALLOCATABLE :: R
!-----------------------------------------------------------------------

lowl = 1.d0                         ! limit of integration low
upl  = 3.d0                         ! limit of integration up

WRITE(*,FMT='(A,E9.3)') 'Using a threshold equal to: ', threshold
WRITE(*,*) '*********************************************'
WRITE(*,*) 'Trapezoidal rule'
CALL Trapezoid(lowl, upl, threshold)

WRITE(*,*) '*********************************************'
WRITE(*,*) 'Simpson'
CALL Simpson(lowl, upl, threshold)

WRITE(*,*) '*********************************************'
WRITE(*,*) 'Romberg'
d = 20
ALLOCATE(R(d,d))
CALL Romberg(R, d, lowl, upl, threshold)

WRITE(*,*) '*********************************************'
WRITE(*,*) 'Gauss-Legendre method'

aux1 = (upl-lowl)/2.d0
aux2 = (upl+lowl)/2.d0
n    = 2
e    = 0.d0; er = 1.d0

DO WHILE ( (er .GT. threshold) .AND. (n .LE. 10) )
  ALLOCATE(t(n), w(n))
  CALL sub_GauLeg(-1.d0, 1.d0, t, w, n)
  s = 0.d0
  DO i = 1, n
    s = s + w(i)*f(aux1*t(i) + aux2)
  ENDDO
  integra = aux1*s
  er      = DABS(integra - e)
  e       = integra
  DEALLOCATE(t, w)             !it will be different dimension
  n = n + 1                    !next iteration
ENDDO

WRITE(*,FMT='(A,F9.7,A,I2)') "The integral value is: ", &
                                       integra, " with quadrature ", n-1

!-----------------------------------------------------------------------
                                CONTAINS
!-----------------------------------------------------------------------
    FUNCTION f(x) RESULT (val)               ! The function to integrate
    IMPLICIT NONE
    REAL*8 :: x, val

    val = DSIN(x*x) - DCOS(2.d0*x)
    ENDFUNCTION f
!-----------------------------------------------------------------------
    SUBROUTINE Trapezoid(lowl,upl,threshold)
    IMPLICIT NONE
    INTEGER :: i, n, iter
    REAL*8  :: lowl, upl, threshold, e, h, p, s, integra, er

    !n how many trapezoids, e and er for threshold
    n = 2; e = 0.d0; er = 1.d0; iter = 0

    WRITE(*,*) "Integral value", "       error"

    DO WHILE (er .GT. threshold)                !until converge
      h = (upl-lowl)/REAL(n)                    !base trapezoid size
      p = (0.5d0*h)*(f(lowl)+f(upl))
      s = 0.d0
      DO i = 1, n-1
        s = s + h*f(lowl+i*h)
      ENDDO
      integra = p+s                      !integral value
      n       = 2*n                      !twice the number of trapezoids
      er      = DABS(e - integra)
      e       = integra
      iter    = iter + 1
      WRITE(*,FMT='(F14.5,F14.5)') integra, er
    ENDDO

    WRITE(*,FMT='((A),F9.7,(A),I4,(A))') "The integral value &
                  is: ", integra, " computed with ", iter, " iterations"

    ENDSUBROUTINE Trapezoid
!-----------------------------------------------------------------------
    SUBROUTINE Simpson(lowl,upl,threshold)
    IMPLICIT NONE
    INTEGER :: i, n
    REAL*8  :: lowl, upl, threshold, e, h, p, s1, s2, integra, er

    !n how many intervals, e and er for threshold
    n = 2; e = 1.d0; er = 1.d0

    WRITE(*,*) "Integral value", "       error"

    DO WHILE (er .GT. threshold)
      h  = (upl - lowl)/n              !base
      p  = f(lowl) + f(upl)
      s1 = 0.d0;  s2 = 0.d0

      DO i = 2, n-2, 2                 !even
        s1 = s1 + f(lowl+i*h)
      ENDDO

      DO i = 1, n-1, 2                 !odd
        s2 = s2 + f(lowl+i*h)
      ENDDO

      integra = (h/3.d0)*(p + 2.d0*s1 + 4.d0*s2)
      n  = 2*n                         !twice the number of trapezoids
      er = DABS(e - integra)
      e  = integra
      WRITE(*,FMT='(F14.5,F14.5)') integra, er
    ENDDO

    WRITE(*,FMT='((A),F9.7,(A),I4,(A))') "The integral value &
                         is: ", integra, " computed with ", n, " points"
    ENDSUBROUTINE Simpson
!-----------------------------------------------------------------------
    SUBROUTINE Romberg(R,d,lowl,upl,threshold)
    IMPLICIT NONE
    INTEGER                 :: i, j, k, n, d
    INTEGER, DIMENSION(2,2) :: flag
    REAL*8                  :: lowl, upl, threshold, e, er, s
    REAL*8                  :: aux1, aux2
    REAL*8, DIMENSION(d)    :: h
    REAL*8, DIMENSION(d,d)  :: R
    CHARACTER(LEN=60)       :: strfor, str

    DO k = 1, d
      h(k) = (upl - lowl)/(2.d0**(k-1))
    ENDDO

    R(1,1) = (h(1)/2.d0)*(f(lowl)+f(upl))

    e = 1.d0; er = 1.d0
    DO k = 2, d          !R_{k,1}   k>1
      s = 0.d0
      DO i = 1, (2**(k-2))
        s = s + f(lowl+(2*i -1)*h(k))
      ENDDO
      R(k,1) = 0.5d0*(R(k-1,1) + h(k-1)*s)
      er = DABS(e - R(k,1))
      e  = R(k,1)
      IF (er .GT. threshold) THEN
        flag(1,1) = k+1
        flag(1,2) = 1
      ENDIF
    ENDDO

    e = 1.d0; er = 1.d0
    j = 2; k = j         !R_{k,j}    k>=j>1

    DO WHILE ( (k .LE. d) .AND. (j .LE. d) )
      aux1   = R(k,j-1) - R(k-1,j-1)
      aux2   = 4.d0**(j-1) - 1.d0
      R(k,j) = R(k,j-1) + aux1/aux2
      er     = DABS(e - R(k,j))
      e      = R(k,j)
      IF (er .GT. threshold) THEN
        flag(2,1) = k+1
        flag(2,2) = j
      ENDIF
      IF (k .EQ. d) THEN
        j = j + 1
        k = j
      ELSE
        k = k + 1
      ENDIF
    ENDDO

    WRITE(*,FMT='(A,F9.7)') 'The integral value is: ', R(flag(1,1), &
                                                              flag(1,2))
    WRITE(*,*) 'And at least, in the positions: '
    DO i = 1, 2
      WRITE(*,FMT='(A3,I2,A1,I2,A1)') 'R_{',flag(i,1),',',flag(i,2),'}'
    ENDDO

    WRITE(*,*)
    WRITE(*,*) 'in: '
    WRITE(*,FMT='(A5,*(A6,I2,A2))') '*****', ('R_{k,',i,'} ', i=1,10)
    DO i = 1, 10
      strfor = '(A2,I2,'//TRIM(str(i))//'(F10.5))'
      WRITE(*,FMT=strfor) 'k=', i, (R(i,j), j=1,i)
    ENDDO
    ENDSUBROUTINE Romberg
!-----------------------------------------------------------------------
    SUBROUTINE sub_GauLeg(X1, X2, t, w, n) !Gauss-Legendreren SUBROUTINE
    ! t and w what gives us value
    ! Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian
    ! Quadrature integration of polynomial functions. For normalized
    ! lower and upper limits of integration: -1.0 & 1.0, and given n,
    ! this routine calculates, arrays t(1:n) and  w(1:n) of length n,
    ! containing the abscissas and weights of the Gauss-Legendre n-point
    ! quadrature formula.
    ! Another way to remove X1 and X2, and XL and XM
    ! And in the end will be t(i) = -z and t(n+1-i)= +z
    IMPLICIT NONE
    INTEGER                           :: m, i, j
    INTEGER, INTENT(IN)               :: n    !Number of Gaussian points
    REAL*8                            :: XM, XL, X1, X2, P1, P2, P3
    REAL*8                            :: Z1, Z, PP
    REAL*8, PARAMETER                 :: EPS = 1.d-14!Relative precision
    REAL*8, PARAMETER                 :: pi = DACOS(-1.d0)  !pi constant
    REAL*8, DIMENSION(n), INTENT(OUT) :: W,t

    m = (n + 1) / 2                 !Roots are symmetric in the interval

    !The coats are going to be X1 = -1 and X2 = 1, Gauss-Legendre 
    XM = 0.5d0*(X1 + X2)
    XL = 0.5d0*(X2 - X1)

    !Loop over the desired roots
    DO i = 1,m
      Z = DCOS (pi * (i - 0.25d0)/(n + 0.5d0))
      !Starting with the above approximation to the i-th root,
      !we enter the main loop of refinement by NEWTON'S method   
 10   P1 = 1.d0
      P2 = 0.d0

      !Loop up the recurrence relation to get the Legendre
      !polynomial evaluated at z                
      DO j = 1,n
        P3 = P2
        P2 = P1
        P1 = ((2.d0 * j - 1.d0) * Z * P2 - (j - 1.d0) * P3)/j
      ENDDO

      !p1 is now the desired Legendre polynomial.
      !We next compute pp, its derivative, by a standard relation
      !involving also p2, the polynomial of one lower order. 
      PP = n * (Z * P1 - P2)/(Z * Z - 1.d0)
      Z1 = Z
      Z  = Z1 - P1/PP                      ! Newton's Method

      IF (DABS(Z-Z1) .GT. EPS) GOTO 10     ! Really old fashion :b

      !Roots will be symmetric about the origin  
      t(i)         = XM - XL * Z
      t(n + 1 - i) = XM + XL * Z
      !Compute the weight and its symmetric counterpart 
      W(i)         = 2.d0 * XL/((1.d0 - Z * Z) * PP * PP)
      W(n + 1 - i) = W(i)
    ENDDO

    RETURN        !not neccesary
    ENDSUBROUTINE Sub_GauLeg
ENDPROGRAM Integration
!***********************************************************************
!---- For fancy style, transform a number to charecter
CHARACTER(LEN=60) FUNCTION str(q)
    INTEGER, INTENT(IN) :: q
    WRITE (str,*) q
    str = ADJUSTL(str)
    ENDFUNCTION str 
