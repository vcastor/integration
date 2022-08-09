PROGRAM monte_carlo

!***********************************************************************
!----     This program computes a integral using the Monte Carlo    ----
!----    method, we will use thepseudo random numbers integraded    ----
!----                        in FORTRAN                             ----
!----                                                               ----
!----                       ~VCastor 2021                           ----
!***********************************************************************
! There is a pdf with documentation about the math theory necessary to
! understand properly the code. MonteCarlo.pdf
!***********************************************************************

IMPLICIT NONE
INTEGER :: i, j, k, l
REAL*8  :: upl, lowl, upy, lowy
!                            f(x) is the function that we will integrate

!***********************************************************************
upl  = 5.d0
lowl = 2.d0
lowy = f(lowl) - 2.d0       !for ordinate-axis
upy  = f(upl)  + 2.d0       !for odrinate-axis

WRITE(*,*) '******************************'
WRITE(*,*) 'Using Sample Mean we get'
CALL SampleMean(lowl,upl)

WRITE(*,*) '******************************'
WRITE(*,*) 'Using Hint and Miss we get'
CALL HintAndMiss(lowl,upl,lowy,upy)


!***********************************************************************
                                  CONTAINS
!***********************************************************************
    FUNCTION f(x) RESULT (value)       ! function that we will integrate
    IMPLICIT NONE
    REAL*8 :: x, val

    val = 2.d0*x*x - 3.d0*x
    ENDFUNCTION f
!-----------------------------------------------------------------------
    SUBROUTINE SampleMean (lowl,upl)
    IMPLICIT NONE
    INTEGER*8                         :: i, j, k, d
    REAL*8                            :: upl, lowl, aux, s, integrate
    REAL*8, DIMENSION(:), ALLOCATABLE :: x
    CHARACTER(LEN=18)                 :: printf
    
    DO k = 1, 8
      d = 10**k
      ALLOCATE(x(d))
      CALL RANDOM_NUMBER(x)      !that random numbers are in range (0,1)
      DO i = 1, d
        !new range (0,{up-low})
        x(i) = x(i)*(upl-lowl) 
        !new range (lowl,upl)
        x(i) = x(i)+lowl
      ENDDO

      s = 0.d0
      DO i = 1, d
        s = s + f(x(i))
      ENDDO
      integrate = (upl-lowl)*s/REAL(d)

      printf = '(A8,I1,A17,F9.6)'
      WRITE(*,FMT=printf) "with 10^", k, " random numbers: ",integrate
      DEALLOCATE(x)
    ENDDO
    ENDSUBROUTINE SampleMean
!-----------------------------------------------------------------------
    SUBROUTINE HintAndMiss(lowl,upl,lowy,upy)
    IMPLICIT NONE
    INTEGER*8                         :: nin, i, j, k, d
    REAL*8                            :: Ar, integrate, lowl, upl
    REAL*8                            ::lowy, upy
    REAL*8, DIMENSION(:), ALLOCATABLE :: x, y, fun
    CHARACTER(LEN=18)                 :: printf

    DO k = 1, 8
      d = 10**k
      ALLOCATE(x(d),y(d),fun(d))
      CALL RANDOM_NUMBER(x)
      CALL RANDOM_NUMBER(y)
      DO i = 1, d
        x(i) = x(i)*(upl-lowl)
        x(i) = x(i)+lowl

        y(i) = y(i)*(upy-lowy)
        y(i) = y(i)+lowy
      ENDDO

      DO i = 1, d
        fun(i) = f(x(i))
      ENDDO

      nin = 0
      DO i = 1, d
        IF ( y(i) .LT. fun(i) ) THEN
          nin = nin + 1
        ENDIF
      ENDDO
        
      Ar        = (upl-lowl)*(upy-lowy)
      integrate = Ar*REAL(nin)/REAL(d)

      printf = '(A8,I1,A17,F9.6)'
      WRITE(*,FMT=printf) "with 10^", k, " random numbers: ",integrate
      DEALLOCATE(x,y,fun)
    ENDDO
    ENDSUBROUTINE HintAndMiss
!-----------------------------------------------------------------------
ENDPROGRAM monte_carlo
