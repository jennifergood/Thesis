MODULE OCEA_MATRIX_LIB
USE PROBLEM_DATA
USE EB_VARIABLE_HANDLING
USE OCEA_MATHLIB
USE OCEA_IO
IMPLICIT NONE

!=======================================================
! LINEAR ALGEBRA ROUTINES SUPPORTED: DOUBLE PRECISION
!=======================================================
! F90_SYM_LINEQN	    CHOLESKEY DRIVER FOR POSITIVE DEFINITE MATRIX
! F90_CHOLESKEY_INV_MATRIX  CHOLESKEY INVERSE FOR SYMMETRIC MATRIX
! F90_CHOLDC                CHOLESKEY FACTORIZATION 
! F90_CHOLSL                CHOLESKEY BACKSUBSTITUTION
! F90_LU_DECOMP_AXeB        LU DECOMPOSITION FOR NON-SYMMETRIC MATRIX
! F90_LU_DECOMP_INV_MATRIX  LU INVERSE FOR NON-SYMMETRIC MATRIX
! F90_ludcmp                LU DECOMPOSITION FACTORIZATION
! F90_LUBKSB                LU DECOMPOSITION BACKSUBSTITUTION
! F90_QRDCMP                QR DECOMPOSITION FACTORIZATION
! F90_QRSOLV                QR DECOMPOSITION BACKSUBSTITUTION
! F90_RSOLV                 QR DECOMPOSITION UTILITIY
! EYE                       DIAGONAL MATRIX  
!====================================================
! LINEAR ALGEBRA ROUTINES SUPPORTED: OCEA EMBEDDED
!====================================================
! EB_CHOLDC                 CHOLESKEY FACTORIZATION
! EB_cholsl                 CHOLESKEY BACKSUBSTITUTION
! EB_CHOL_DECOMP_AXeB       CHOLESKEY LINEAR EQN
! EB_CHOL_DECOMP_INV_MAT    CHOLESKEY MATRIX INVERSE
! EB_LUDCMP                 LU FACTORIZATION
! EB_LUBKSB                 LU BACKSUBSTITUTION
! EB_LU_DECOMP_AXeB         LU LINEAR EQN
! EB_LU_DECOMP_INV_MATRIX   LU MATRIX INVERSE
! EB_AXeB                   GAUSSIAN ELIMINATION NO PIVITING LIN EQN
! EB_MAT_INV_GELIM          GAUSSIAN ELIMINATION NO PIVITING MAT INV
! EB_JACOBI                 JACOBI ROTATION FOR SYM EIGENVALUE/VECTOR



CONTAINS

!=======================================================================
!     CHOLESKEY DECOMPOSITION
!=======================================================================
      SUBROUTINE f90_SYM_LINEQN( A, B, XSOLN)
!	  THIS PROGRAM SOLVES THE LINEAR  MATRIX EQUATION AX = B FOR X.  THE
!	  INPUT MATRIX IS ASSUMED SYMMETRIC AND POSITIVE DEFINITE.  THE SOLUTION
!	  IS OBTAINED BY  USING CHOLESKEY DECOMPOSITION.
!=================================================
! COPYRIGHT(C) 2003 JAMES D. TURNER, AMDYN SYSTEMS
!=================================================	
      IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES

      REAL(DP), DIMENSION(:,:),       INTENT(IN ):: A
      REAL(DP), DIMENSION(SIZE(A,2)), INTENT(IN ):: B
      REAL(DP), DIMENSION(SIZE(A,2)), INTENT(OUT):: XSOLN

!.....LOCAL VARIABLES

      REAL(DP), DIMENSION(SIZE(A,2),SIZE(A,2)):: ASOLVE
      REAL(DP), DIMENSION(SIZE(A,2)          ):: P
      REAL(DP), DIMENSION(SIZE(A,2)          ):: ERROR
      INTEGER:: I,J
      LOGICAL:: PRT

      COPY_A_MATRIX: IF( SIZE(A,1) /= SIZE(A,2) ) THEN
          DO I = 1, SIZE(A,2)
             ASOLVE(I,:) = A(I,:)
          END DO
      ELSE
          ASOLVE = A
      END IF COPY_A_MATRIX

      CALL F90_CHOLDC( ASOLVE, P)

      CALL F90_CHOLSL( ASOLVE, P, B, XSOLN)

      PRT=.FALSE.
      PRINT_CK_SOLN_ACCURACY: IF( PRT ) THEN
          ERROR = MATMUL( A,XSOLN) - B
          WRITE(6,*)
          WRITE(6,*)' MODULE:F90_LINEAR_ALGEBRA--ROUTINE F90_SYM_LINEQN'
          WRITE(6,*)' CHECKING SOLUTION ACCURACY'
          WRITE(6,*)
          WRITE(6,'( A/,1000(4E12.3/))')' ERROR X', &
                   ( ERROR(I), XSOLN(I), I=1,SIZE(A,2))
          WRITE(6,*)
          WRITE(6,*)' FINISHED IN SYM_LINEQN'
      END IF PRINT_CK_SOLN_ACCURACY

      END subroutine f90_sym_lineqn

!=========================================================================
SUBROUTINE F90_CHOLESKEY_INV_MATRIX( A, INV_A )
!	THIS PROGRAM COMPUTES THE INVERSE OF A SYM. MATRIX USING CHOLESKEY 
!	DECOMPOSITION.  THE MATRIX IS ASSUMED TO BE
!	SYMMETRIC AND POSITIVE-DEFINITE.
!
!.....INPUT
!	A		SYMMETRIC MATRIX
!
!.....OUTPUT
!	INV_A	INVERSE OF THE SYMMETRIC INPUT MATRIX
!
!.....EXTERNAL PROCEDURES
!	F90_CHOLDC	VERSION OF CHOLESKEY DECOMPOSITION 
!	F90_CHOLSL	VERSION OF CHOLESKEY BACKSUBSTITUTION
!====================================
! COPYRIGHT(C) 2003 JAMES D. TURNER
!====================================
USE PROBLEM_DATA
IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES
REAL(DP),DIMENSION(:,:),INTENT(IN ):: A
REAL(DP),DIMENSION(:,:),INTENT(OUT):: INV_A 
!.....LOCAL VARIABLES
REAL(DP),DIMENSION(SIZE(A,2),SIZE(A,2)):: ASAVE
REAL(DP),DIMENSION(SIZE(A,2)          ):: P
REAL(DP),DIMENSION(SIZE(A,2)          ):: COL
INTEGER:: I,N
N=SIZE(A,1)
COL = 0.0D0; COL(1) = 1.0D0                             ! DEFINE COL OF IDENTITY MATRIX

COPY_A_MATRIX:IF( SIZE(A,1) /= SIZE(A,2) ) THEN
    DO I = 1, SIZE(A,2)
       ASAVE(I,:) = A(i,:)
    END DO
ELSE
    ASAVE = A
END IF COPY_A_MATRIX                              ! SAVE INITIAL INPUT MATRIX

CALL F90_CHOLDC( ASAVE, P        )

CALL F90_CHOLSL( ASAVE, P, COL, INV_A(:,1) )      ! COMPUTE FIRST COLUMN OF INVERSE MATRIX

DO I = 2, SIZE(A,2)

   COL(I-1) = 0.0D0; COL(I  ) = 1.0D0                   ! ADVANCE NEXT COL OF IDENTITY MATRIX

   CALL F90_CHOLSL( ASAVE, P, COL, INV_A(:,I) )   ! COMPUTE ITH COLUMN OF INVERSE MATRIX

END DO

END SUBROUTINE F90_CHOLESKEY_INV_MATRIX

!=======================================================================
!   THIS PROGRAM COMPUTES A CHOLESKY DECOMPOSITION OF A SYMMETRIX MATRIX
!
SUBROUTINE F90_CHOLDC( a, p)
      IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES

      REAL(DP), DIMENSION( :, :), INTENT(INOUT):: A
      REAL(DP), DIMENSION( :   ), INTENT(OUT  ):: P

!.....LOCAL VARIABLES

      REAL(DP):: SUM
      INTEGER:: i,j

      do i = 1,SIZE(A,2)
        do j = i,SIZE(A,2)
          sum = a(i,j) - dot_product( a(i,i-1:1:-1), a(j,i-1:1:-1) )
          if(i==j)then

             if(sum <= 0.0D0)then
                write(6,*)
                write(6,*)'===================================================='
                write(6,*)'MODULE:F90_LINEAR_ALGEBRA--ROUTINE F90_choldc failed'
                write(6,*)'===================================================='
               !====================================================
               ! don't print a matrix here: lower part transformed!!
               !====================================================
                stop'choldc input matrix not positive definite'
             end if

             p(i)=Dsqrt(sum)
          else
             a(j,i)=sum/p(i)
          end if

        END DO !J
      END DO !I

      END SUBROUTINE F90_CHOLDC
!  (C) Copr. 1986-92 Numerical Recipes Software

!=======================================================================
!   THIS PROGRAM COMPUTES A BACK SUBSTITUTION GIVEN A HAS BEEN DECOMPOSED
!   BY CALLING CHOLDC TO PROVIDE SOLUTION FOR A SYMMETRIC LINEAR EQUATION

SUBROUTINE F90_CHOLSL( a, p, b, x)
      USE PROBLEM_DATA
      IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES

      REAL(DP), DIMENSION(:,:), INTENT(IN ):: A
      REAL(DP), DIMENSION(:  ), INTENT(IN ):: P, B
      REAL(DP), DIMENSION(:  ), INTENT(OUT):: X

!.....LOCAL VARIABLES

      INTEGER::  i,N
      REAL(DP):: Sum,tmp
      N=SIZE(A,2)

      do i=1,N                  ! Solve L.y = b, storing y in x
         x(i) = ( b(i) - dot_product( a(i,i-1:1:-1), x(i-1:1:-1) ) )/p(i)
      end do
      do i=N,1,-1               ! Solve LT.x = y
         x(i) = ( x(i) - dot_product( a(i+1:n,i), x(i+1:n) ) )/p(i)
      end do

      END SUBROUTINE F90_CHOLSL
!  (C) Copr. 1986-92 Numerical Recipes Software

!===============================================
!     LU DECOMPOSITION
!===============================================
SUBROUTINE F90_LU_DECOMP_AXeB( A, B, X )
    USE PROBLEM_DATA
    uSE OCEA_IO
    IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES

    REAL(DP),DIMENSION(:,:),INTENT(IN ):: A
    REAL(DP),DIMENSION(:  ),INTENT(IN ):: B
    REAL(DP),DIMENSION(:  ),INTENT(OUT):: X

!.....LOCAL VARIABLES

    REAL(DP),DIMENSION(SIZE(A,2),SIZE(A,2)):: ASAVE 
    INTEGER, DIMENSION(SIZE(A,2))          :: INDX  
    REAL(DP),DIMENSION(SIZE(A,2))          :: BSAVE,TMP
    REAL(DP):: D
    INTEGER:: I,J
    LOGICAL:: PRT

    X = B
    COPY_A_MATRIX:IF ( SIZE(A,1) /= SIZE(A,2) )THEN
       DO I = 1, SIZE(A,2)
          ASAVE(I,:) = A(I,:)
       END DO !I
    ELSE
       ASAVE = A
    END IF COPY_A_MATRIX

    INDX = 0

    CALL F90_LUDCMP( ASAVE, INDX, D )

    CALL F90_LUBKSB( ASAVE, INDX, X )

    PRT=.FALSE.
    IF(PRT) THEN

       TMP=MATMUL( A, X ) - B

       WRITE(6,*)'F90_LU_DECOMP__AXeB:ERROR CHECK'
       CALL OUT_MATRIX( 6, TMP, SIZE(A,2), SIZE(A,2), 1)
    END IF

END SUBROUTINE F90_LU_DECOMP_AXeB

!====================================================
SUBROUTINE f90_LU_DECOMP_INV_MATRIX( A, INV_A )
!	THIS PROGRAM COMPUTES THE INVERSE OF A GEN. MATRIX USING LU 
!	DECOMPOSITION.  THE MATRIX IS ASSUMED TO BE
!	GENERAL AND NON-POSITIVE-DEFINITE.
!
!.....INPUT
!	A		GENERAL MATRIX
!
!.....OUTPUT
!	INV_A	INVERSE OF THE GENERAL INPUT MATRIX
!
!.....EXTERNAL PROCEDURES
!	f90_LUDCMP	VERSION OF LU DECOMPOSITION 
!	f90_LUBKSB	VERSION OF LU BACKSUBSTITUTION
!====================================
! COPYRIGHT(C) 2003 JAMES D. TURNER
!====================================    
    USE PROBLEM_DATA
    IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES

    REAL(DP),DIMENSION(:,:),INTENT(IN ):: A
    REAL(DP),DIMENSION(:,:),INTENT(OUT):: INV_A 

!.....LOCAL VARIABLES

    REAL(DP),DIMENSION(SIZE(A,2),SIZE(A,2)):: ASAVE
    INTEGER :: N,I
    REAL(DP):: D
    INTEGER,DIMENSION(SIZE(A,1)):: INDX

    COPY_A_MATRIX:IF ( SIZE(A,1) /= SIZE(A,2) )THEN
       DO I = 1, SIZE(A,2)
          ASAVE(I,:) = A(I,:)
       END DO !I
    ELSE
       ASAVE = A
    END IF COPY_A_MATRIX

    N = SIZE(A,2)

    INV_A = 0.0D0

    DO I = 1, N
       INV_A(I,I) = 1.0D0
    END DO

    CALL f90_LUDCMP( ASAVE, INDX, D )

    DO I = 1, N
       CALL f90_LUBKSB( ASAVE, INDX, INV_A(:,I) )
    END DO

END SUBROUTINE f90_LU_DECOMP_INV_MATRIX

!==========================================================
SUBROUTINE F90_ludcmp( a, indx ,d )
    USE PROBLEM_DATA
    IMPLICIT NONE

    REAL(DP),PARAMETER                     :: TINY=1.0e-20
    INTEGER, PARAMETER                     :: NMAX = 500

!.....ARGUMENT LIST VARIABLES

    REAL(DP),DIMENSION(:,:  ),INTENT(INOUT):: A
    INTEGER, DIMENSION(:    ),INTENT(INOUT):: indx

!.....LOCAL VARIABLES

    REAL(DP),                 INTENT(OUT  ):: d
    INTEGER                                :: i,imax,j,k
    REAL(DP)                               :: aamax,dum,sum
    REAL(DP),DIMENSION( SIZE(A,1) )        :: vv, ROW
    REAL(DP)                               :: TMP
    INTEGER                                :: N

    N=SIZE(A,2)

      d=ONE

      VV = MAXVAL( ABS(A), DIM = 2 )
      IF ( ANY(VV==ZERO )) STOP "SINGULAR MATRIX IN LUDCMP"
      VV = ONE/VV

      do 19 j=1,n

        do 14 i=1,j-1
          a(i,j)= A(I,J) - DOT_PRODUCT( A(I, 1:I-1), A(1:I-1, J) )
14      continue

        aamax=ZERO
        do 16 i=j,n
          a(i,j)= A(I,J) - DOT_PRODUCT( A(I, 1:J-1), A(1:J-1, J) )
          dum=vv(i)*abs( A(I,J) )
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue

        if (j.ne.imax)then

          ROW = A(IMAX,:); A(IMAX,:) = A(J,:); A(J,:) = ROW

          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.ZERO)a(j,j)=TINY
        if(j.ne.n)then

          dum=ONE/a(j,j)
          A(J+1:N, J) = A(J+1:N, J)*DUM

        endif
19    continue

      END SUBROUTINE F90_LUDCMP
!  (C) Copr. 1986-92 Numerical Recipes Software 

!=====================================================
SUBROUTINE F90_lubksb(a,indx,b)
      USE PROBLEM_DATA
      IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES

      REAL(DP),DIMENSION(:,:),INTENT(INOUT):: A
      INTEGER, DIMENSION(:  ),INTENT(INOUT):: indx
      REAL(DP),DIMENSION(:  ),INTENT(INOUT):: B

!.....LOCAL VARIABLES

      INTEGER:: i,ii,j,ll
      REAL(DP):: sum
      INTEGER::  N
      N=SIZE(A,2)
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          SUM = SUM - DOT_PRODUCT( A(I,II:I-1), B(II:I-1) )
        else if (sum.ne.ZERO) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
         b(i) = ( B(I) - DOT_PRODUCT( A(I,I+1:N), B(I+1:N) ) )/A(I,I)
14    continue
      
END SUBROUTINE F90_LUBKSB
!  (C) Copr. 1986-92 Numerical Recipes Software 

!=======================================
!     QR DECOMPOSITION
!=======================================
SUBROUTINE F90_QRDCMP(a,n,np,c,d,sing)
      USE PROBLEM_DATA
      IMPLICIT NONE

      INTEGER n,np
      REAL(DP) a(np,np),c(n),d(n)
      LOGICAL sing
      INTEGER i,j,k
      REAL(DP) scale,sigma,sum,tau
      sing = .false.
      scale = ZERO
      do 17 k = 1,n-1
        do 11 i = k,n
          scale = max(scale,abs(a(i,k)))
11      continue
        if(scale.eq.ZERO)then
          sing = .true.
          c(k) = ZERO
          d(k) = ZERO
        else
          do 12 i = k,n
            a(i,k) = a(i,k)/scale
12        continue
          sum = ZERO
          do 13 i = k,n
            sum = sum + a(i,k)**2
13        continue
          sigma = sign(sqrt(sum),a(k,k))
          a(k,k) = a(k,k) + sigma
          c(k) = sigma*a(k,k)
          d(k) = -scale*sigma
          do 16 j = k + 1,n
            sum = ZERO
            do 14 i = k,n
              sum = sum + a(i,k)*a(i,j)
14          continue
            tau = sum/c(k)
            do 15 i = k,n
              a(i,j) = a(i,j)-tau*a(i,k)
15          continue
16        continue
        endif
17    continue
      d(n) = a(n,n)
      if(d(n).eq.ZERO)sing = .true.
  
END SUBROUTINE F90_QRDCMP
!  (C) Copr. 1986-92 Numerical Recipes Software 

!============================================      
SUBROUTINE F90_qrsolv(a,n,np,c,d,b)
      USE PROBLEM_DATA
      IMPLICIT NONE

      INTEGER n,np
      REAL(DP) a(np,np),b(n),c(n),d(n)
!.....USES rsolv
      INTEGER i,j
      REAL(DP) sum,tau
      do 13 j=1,n-1
        sum=0.0D0
        do 11 i=j,n
          sum=sum+a(i,j)*b(i)
11      continue
        tau=sum/c(j)
        do 12 i=j,n
          b(i)=b(i)-tau*a(i,j)
12      continue
13    continue
      call F90_rsolv(a,n,np,d,b)
      return
      END SUBROUTINE F90_QRSOLV
!  (C) Copr. 1986-92 Numerical Recipes Software

!==============================================
SUBROUTINE F90_rsolv(a,n,np,d,b)
USE PROBLEM_DATA
IMPLICIT NONE

      INTEGER,                  INTENT(IN   ):: n,np
      REAL(DP),DIMENSION(NP,NP),INTENT(IN   ):: a
      REAL(DP),DIMENSION(N    ),INTENT(IN   ):: d
      REAL(DP),DIMENSION(N    ),INTENT(INOUT):: b
      INTEGER:: i
      b(n)=b(n)/d(n)
      do i=n-1,1,-1
         b(i)=(b(i)-DOT_PRODUCT( A(I,I+1:N), B(I+1:N) ) )/d(i)
      END DO !I
END SUBROUTINE F90_RSOLV
!  (C) Copr. 1986-92 Numerical Recipes Software 


!=============================================
FUNCTION EYE(DIM,VALUE)
    USE PROBLEM_DATA
    IMPLICIT NONE

    INTEGER, INTENT(IN):: DIM
    REAL(DP), DIMENSION(DIM,DIM):: EYE
    REAL(DP),INTENT(IN):: VALUE
    INTEGER:: I
    EYE   = 0.0D0
    DO I=1,DIM
       EYE(I,I) = VALUE
    END DO
END FUNCTION EYE


!============================================
!     OCEA EMBEDDED DOUBLE PRECISION
!============================================

!============================================
SUBROUTINE EB_choldc( a, p )
!	THIS CODE IS A MODIFIED VERSION OF NUMERICAL RECIPES SOFTWARE
!	PROGRAM CHOLDC FOR COMPUTING A CHOLESKEY DECOMPOSITION.  THE
!	MODIFIED VERSION NOW SUPPORTS OCEA ALGEBRA FOR GENERATING
!	FIRST AND SECOND ORDER PARTIAL DERIVATIVE CALCULAIONS FOR
!	THE CHOLESKEY DECOMPOSITION OF AN OCEA SYMMETRIC MATRIX
!
!.....INPUT
!
!	A	OCEA SYMMETRIC MATRIX
!
!		= [ A, DEL(A), DEL(DEL(A)) ].
!
!.....OUTPUT
!	P	OCEA VECTOR OF DIAGONAL ELEMENT FOR CHOLESKY FACTOR L.
!
!		= [ P, DEL(P), DEL(DEL(B)) ].
!
!===============================================================
!	MODIFICATION AUTHOR: JAMES D. TURNER, AMDYN SYSTEMS, 2003
!	CONTACT:             JAMESDANIELTURNER@HOTMAIL.COM
!===============================================================
      USE EB_VARIABLE_HANDLING
      IMPLICIT NONE

!...ARGUMENT LIST VARIABLES

      TYPE(EB), DIMENSION(:,:), INTENT(INOUT):: A
      TYPE(EB), DIMENSION(:  ), INTENT(OUT  ):: P

!...LOCAL VARIABLES

      INTEGER:: n, i, j, k
      TYPE(EB):: sum

      N = SIZE(A,2)
      do 13 i=1,n
        do 12 j=i,n
          sum=a(i,j)
          do 11 k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
11        continue
          if(i.eq.j)then
            if(sum%E.le.0.0D0)STOP 'EB_choldc failed'
            p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
12      continue
13    continue
    
END SUBROUTINE EB_choldc
!  (C) Copr. 1986-92 Numerical Recipes Software

!==============================================
SUBROUTINE EB_cholsl( a, p, b, x )
!	THIS CODE IS A MODIFIED VERSION OF NUMERICAL RECIPES SOFTWARE
!	PROGRAM CHOLSL FOR COMPUTING A CHOLESKEY DECOMPOSITION BACKSUBSTITUTION.
!	THE MODIFIED VERSION NOW SUPPORTS OCEA ALGEBRA FOR GENERATING
!	FIRST AND SECOND ORDER PARTIAL DERIVATIVE CALCULAIONS FOR
!	THE CHOLESKEY DECOMPOSITION OF AN OCEA SYMMETRIC MATRIX
!
!.....INPUT
!
!	A	OCEA SYMMETRIC MATRIX
!
!		= [ A, DEL(A), DEL(DEL(A)) ].
!
!	B	OCEA RIGHT HAND SIDE (AX = B)
!
!		= [ B, DEL(B), DEL(DEL(B)) ]
!
!.....OUTPUT
!	X	OCEA SOLUTION VECTOR FOR AX = B
!
!		= [ X, DEL(X), DEL(DEL(X)) ]
!
!===============================================================
!	MODIFICATION AUTHOR: JAMES D. TURNER, AMDYN SYSTEMS, 2003
!	CONTACT:             JAMESDANIELTURNER@HOTMAIL.COM
!===============================================================
USE EB_VARIABLE_HANDLING
IMPLICIT NONE

!...ARGUMENT LIST VARIABLES
    TYPE(EB), DIMENSION(:,:),     INTENT(IN   ):: A
    TYPE(EB), DIMENSION(:  ),     INTENT(IN   ):: B, P
    TYPE(EB), DIMENSION(SIZE(B)), INTENT(INOUT):: X
!...LOCAL VARIABLES
    INTEGER:: n, i, k
    INTEGER:: DIM
    TYPE(EB):: sum

   N = DIM
   do 12 i=1,n
      sum=b(i)
      do 11 k=i-1,1,-1
         sum=sum-a(i,k)*x(k)
11       continue
      x(i)=sum/p(i)
12    continue
   do 14 i=n,1,-1
      sum=x(i)
      do 13 k=i+1,n
         sum=sum-a(k,i)*x(k)
13       continue
      x(i)=sum/p(i)
14    continue

END SUBROUTINE EB_cholsl
!  (C) Copr. 1986-92 Numerical Recipes Software 

!==============================================
SUBROUTINE EB_CHOL_DECOMP_AXeB( A, B, X )
!	THIS PROGRAM SOLVES THE LINEAR MATRIX EQUATION OF AN OCEA MATRIX 
!	AND OCEA RHS USING CHOLESKEY DECOMPOSITION.  THE MATRIX%E PART OF 
!	THE MATRIX IS ASSUMED TO BE SYMMETRIC AND POSITIVE-DEFINITE.
!
!.....INPUT
!	A	MATRIX	SYMMETRIC EMBEDDED MATRIX
!
!		= [     M,       DEL(M),      DEL(DEL(M)) ]
!
!	B	RIGHT HAND SIDE OF LINEAR EQUATION
!
!		= [     B,       DEL(B),      DEL(DEL(B)) ]
!
!.....OUTPUT
!	X	LINEAR EQUATION SOLUTION EMBEDDED PROBLEM INPUTS
!
!		= [ INV(M)B, DEL(INV(M)B), DEL(DEL(INV(M)B)) ]
!
!.....EXTERNAL PROCEDURES
!		EB_CHOLDC	EMBEDDED VERSION OF CHOLESKEY DECOMPOSITION 
!		EB_CHOLSL	EMBEDDED VERSION OF CHOLESKEY BACKSUBSTITUTION
!====================================
! COPYRIGHT(C) 2003 JAMES D. TURNER
!====================================	
USE PROBLEM_DATA
USE EB_VARIABLE_HANDLING
IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES

TYPE(EB), DIMENSION(:,:),     INTENT(IN ):: A
TYPE(EB), DIMENSION(:  ),     INTENT(IN ):: B
TYPE(EB), DIMENSION(size(b)), INTENT(OUT):: X

!.....LOCAL VARIABLES

TYPE(EB), DIMENSION(size(b),size(b)):: ATMP 
TYPE(EB), DIMENSION(size(b))        :: P  
TYPE(EB), DIMENSION(size(b))        :: BTMP
TYPE(EB)                            :: D
INTEGER:: I

BTMP = B

COPY_A_MATRIX:IF ( SIZE(A,1) /= SIZE(A,2) )THEN
    DO I = 1, SIZE(A,2)
        ATMP(I,:) = A(I,:)
    END DO !I
ELSE
    ATMP = A
END IF COPY_A_MATRIX

!.....chol-DECOMP FOR LINEAR EQN. SOLN

      CALL EB_CHOLDC( ATMP, P        )

      CALL EB_CHOLSL( ATMP, P, BTMP, X )

END SUBROUTINE EB_CHOL_DECOMP_AXeB

!===============================================
SUBROUTINE EB_CHOL_DECOMP_INV_MAT( A, INV_A )
!	THIS PROGRAM COMPUTES THE INVERSE OF AN OCEA MATRIX USING CHOLESKEY 
!	DECOMPOSITION.  THE MATRIX%E PART OF THE MATRIX IS ASSUMED TO BE
!	SYMMETRIC AND POSITIVE-DEFINITE.
!
!.....INPUT
!	A		SYMMETRIC EMBEDDED MATRIX
!
!			= [     M,       DEL(M),      DEL(DEL(M)) ]
!.....OUTPUT
!	INV_A    	INVERSE OF THE SYMMETRIC EMBEDDED INPUT MATRIX
!
!			= [ INV(M), DEL(INV(M)), DEL(DEL(INV(M))) ]
!
!.....EXTERNAL PROCEDURES
!		EB_CHOLDC	EMBEDDED VERSION OF CHOLESKEY DECOMPOSITION 
!		EB_CHOLSL	EMBEDDED VERSION OF CHOLESKEY BACKSUBSTITUTION
!====================================
! COPYRIGHT(C) 2003 JAMES D. TURNER
!====================================
USE EB_VARIABLE_HANDLING
IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES

TYPE(EB),DIMENSION(:,:),INTENT(IN ):: A
TYPE(EB),DIMENSION(:,:),INTENT(OUT):: INV_A 

!.....LOCAL VARIABLES

TYPE(EB),DIMENSION(SIZE(A,2),SIZE(A,2)):: ASAVE,TMP
TYPE(EB),DIMENSION(SIZE(A,2))          :: P
TYPE(EB),DIMENSION(SIZE(A,2))          :: COL
INTEGER:: I, DIM

DIM = SIZE(A,2)

COPY_A_MATRIX:IF ( SIZE(A,1) /= SIZE(A,2) )THEN
   DO I = 1, DIM
      ASAVE(I,:) = A(I,:)
   END DO !I
ELSE
   ASAVE = A
END IF COPY_A_MATRIX

COL = .constant.0.0D0; COL(1)%E = 1.0D0                      ! DEFINE COL OF IDENTITY MATRIX

CALL EB_CHOLDC( ASAVE, P        )                  ! CHOLESKEY FACTOR INPUT MATRIX

CALL EB_CHOLSL( ASAVE, P, COL, TMP(:,1) )        ! COMPUTE FIRST COL OF INV MATRIX

DO I = 2, DIM

   COL(I-1)%E = 0.0D0; COL(I  )%E = 1.0D0          ! ADVANCE NEXT COL OF IDENT MATRIX

   CALL EB_CHOLSL( ASAVE, P, COL, TMP(:,I) )     ! COMPUTE ITH COLUMN OF INV MATRIX

END DO

COPY_INV_A_MATRIX:IF ( SIZE(A,1) /= SIZE(A,2) )THEN
   INV_A(1:DIM,1:DIM) = TMP(1:DIM,1:DIM)
END IF COPY_INV_A_MATRIX

END SUBROUTINE EB_CHOL_DECOMP_INV_MAT

!==============================================
SUBROUTINE EB_ludcmp(a,indx,d)
!	THIS CODE IS A MODIFIED VERSION OF NUMERICAL RECIPES SOFTWARE
!	PROGRAM CHOLDC FOR COMPUTING A CHOLESKEY DECOMPOSITION.  THE
!	MODIFIED VERSION NOW SUPPORTS OCEA ALGEBRA FOR GENERATING
!	FIRST AND SECOND ORDER PARTIAL DERIVATIVE CALCULAIONS FOR
!	THE CHOLESKEY DECOMPOSITION OF AN OCEA SYMMETRIC MATRIX
!
!.....INPUT
!		A	OCEA SYMMETRIC MATRIX
!
!			= [ A, DEL(A), DEL(DEL(A)) ].
!
!.....OUTPUT
!		P	OCEA VECTOR OF DIAGONAL ELEMENT FOR CHOLESKY FACTOR L.
!
!			= [ P, DEL(P), DEL(DEL(B)) ].
!
!===============================================================
!	MODIFICATION AUTHOR: JAMES D. TURNER, AMDYN SYSTEMS, 2003
!	CONTACT:             JAMESDANIELTURNER@HOTMAIL.COM
!===============================================================
USE EB_VARIABLE_HANDLING
IMPLICIT NONE

!...PARAMETER VARIABLES
REAL(DP),PARAMETER                       :: TINY = 1.0e-20
INTEGER, PARAMETER                       :: NMAX = 500
!...ARGUMENT LIST VARIABLES
TYPE(EB),DIMENSION(:,:),      INTENT(INOUT):: A
INTEGER, DIMENSION(SIZE(A,2)),INTENT(OUT  ):: indx
REAL(DP),                     INTENT(INOUT):: d
!...LOCAL VARIABLES
INTEGER                                  :: N,i,imax,j,k
TYPE(EB)                                 :: aamax,dum,sum
TYPE(EB),DIMENSION(NMAX)                 :: vv
TYPE(EB),DIMENSION(SIZE(A,2))            :: ROW
INTEGER                                  :: DIM

      N = SIZE(A,2)
      d = ONE
      do 12 i = 1,n
        aamax = .constant.ZERO
        do 11 j = 1,n
          if (abs( a(i,j)%E ).gt.aamax%E) aamax = abs( a(i,j) )
11      continue
        if (aamax%E.eq.ZERO) STOP 'singular matrix in ludcmp'
        vv(i) = ONE/aamax
12    continue
      do 19 j = 1,n
        do 14 i = 1,j-1
          sum = a(i,j)
          do 13 k = 1,i-1
            sum = sum-a(i,k)*a(k,j)
13        continue
          a(i,j) = sum
14      continue
        aamax = .constant.ZERO
        do 16 i = j,n
          sum = a(i,j)
          do 15 k = 1,j-1
            sum = sum-a(i,k)*a(k,j)
15        continue
          a(i,j) = sum
          dum = vv(i)*abs( sum )
          if (dum%E.ge.aamax%E) then
            imax = i
            aamax = dum
          endif
16      continue
        if (j.ne.imax)then
            ROW = A(IMAX,:); A(IMAX,:) = A(J,:); A(J,:) = ROW
            d = -d
            vv(imax) = vv(j)
        endif
        indx(j) = imax
        if(a(j,j)%E.eq.ZERO)a(j,j)%E = TINY   
        if(j.ne.n)then
          dum = ONE/a(j,j)
          do 18 i = j+1,n
            a(i,j) = a(i,j)*dum
18        continue
        endif
19    continue

END SUBROUTINE EB_LUDCMP
!  (C) Copr. 1986-92 Numerical Recipes Software 


!==============================================
SUBROUTINE EB_lubksb(a,indx,b)
!	THIS CODE IS A MODIFIED VERSION OF NUMERICAL RECIPES SOFTWARE
!	PROGRAM LUBKDB FOR COMPUTING A LU DECOMPOSITION BACKSUBSTITUTION.
!	THE MODIFIED VERSION NOW SUPPORTS OCEA ALGEBRA FOR GENERATING
!	FIRST AND SECOND ORDER PARTIAL DERIVATIVE CALCULAIONS FOR
!	THE CHOLESKEY DECOMPOSITION OF AN OCEA SYMMETRIC MATRIX
!
!.....INPUT
!	A	OCEA SYMMETRIC MATRIX
!
!		= [ A, DEL(A), DEL(DEL(A)) ].
!
!	B	OCEA RIGHT HAND SIDE (AX = B)
!
!		= [ B, DEL(B), DEL(DEL(B)) ]
!
!.....OUTPUT
!	X	OCEA SOLUTION VECTOR FOR AX = B  (B OVER WRITTEN WITH X ON OUTPUT)
!
!		= [ X, DEL(X), DEL(DEL(X)) ].
!
!===============================================================
!	MODIFICATION AUTHOR: JAMES D. TURNER, AMDYN SYSTEMS, 2003
!	CONTACT:             JAMESDANIELTURNER@HOTMAIL.COM
!===============================================================
      USE EB_VARIABLE_HANDLING
      IMPLICIT NONE

      TYPE(EB), DIMENSION(:,:),      INTENT(IN   ):: a
      INTEGER,  DIMENSION(SIZE(A,2)),INTENT(IN   ):: indx
      TYPE(EB), DIMENSION(SIZE(A,2)),INTENT(INOUT):: b
      INTEGER :: N,i,ii,j,ll
      TYPE(EB):: sum

      N = SIZE(A,2)
      ii = 0
      do 12 i = 1,n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if (ii.ne.0)then
          do 11 j = ii,i-1
            sum = sum-a(i,j)*b(j)
11        continue
        else if (sum%E.ne.ZERO) then
          ii = i
        endif
        b(i) = sum
12    continue
      do 14 i = n,1,-1
        sum = b(i)
        do 13 j = i+1,n
          sum = sum-a(i,j)*b(j)
13      continue
        b(i) = sum/a(i,i)
14    continue

END SUBROUTINE EB_lubksb
!  (C) Copr. 1986-92 Numerical Recipes Software 


!==============================================
SUBROUTINE EB_LU_DECOMP_AXeB( A, B, X )
!	THIS PROGRAM SOLVES THE LINEAR MATRIX EQUATION OF AN OCEA MATRIX 
!	AND OCEA RHS USING LU DECOMPOSITION.  THE MATRIX%E PART OF 
!	THE MATRIX IS ASSUMED TO BE SQUARE AND NON-SYMMETRIC.
!
!.....INPUT
!	MATRIX		NON - SYMMETRIC EMBEDDED MATRIX
!
!			= [     M,       DEL(M),      DEL(DEL(M)) ]
!
!	B		RIGHT HAND SIDE OF LINEAR EQUATION
!
!			= [     B,       DEL(B),      DEL(DEL(B)) ]
!
!
!.....OUTPUT
!	X		LINEAR EQUATION SOLUTION EMBEDDED PROBLEM INPUTS
!
!			= [ INV(M)B, DEL(INV(M)B), DEL(DEL(INV(M)B)) ]
!
!.....EXTERNAL PROCEDURES
!	EB_LUDCMP	EMBEDDED VERSION OF CHOLESKEY DECOMPOSITION 
!	EB_LUBKSB	EMBEDDED VERSION OF CHOLESKEY BACKSUBSTITUTION
!=================================================
! COPYRIGHT(C) 2003 JAMES D. TURNER, AMDYN SYSTEMS
!=================================================	
USE PROBLEM_DATA
USE EB_VARIABLE_HANDLING
IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES

TYPE(EB),DIMENSION(:,:),    INTENT(IN ):: A
TYPE(EB),DIMENSION(:),      INTENT(IN ):: B
TYPE(EB),DIMENSION(SIZE(B)),INTENT(OUT):: X

!.....LOCAL VARIABLES

TYPE(EB),DIMENSION(SIZE(B),SIZE(B)):: ATMP 
INTEGER, DIMENSION(SIZE(B))        :: INDX  
REAL(DP):: D
INTEGER:: I

X = B; INDX  =  0

COPY_A_MATRIX:IF ( SIZE(A,1) /= SIZE(A,2) )THEN
   DO I = 1, SIZE(A,2)
      ATMP(I,:) = A(I,:)
   END DO !I
ELSE
   ATMP = A
END IF COPY_A_MATRIX

!.....LU-DECOMP FOR LINEAR EQN. SOLN

    CALL EB_LUDCMP( ATMP, INDX, D )

    CALL EB_LUBKSB( ATMP, INDX, X )

END SUBROUTINE EB_LU_DECOMP_AXeB

!==============================================
SUBROUTINE EB_LU_DECOMP_INV_MATRIX( A, INV_A )
!	THIS PROGRAM COMPUTES THE INVERSE OF AN OCEA MATRIX USING LU 
!	DECOMPOSITION.  THE MATRIX%E PART OF THE MATRIX IS ASSUMED TO BE
!	NON - SYMMETRIC.
!
!.....INPUT
!	A		NON - SYMMETRIC EMBEDDED MATRIX
!
!			= [     M,       DEL(M),      DEL(DEL(M)) ]
!.....OUTPUT
!	INV_A	INVERSE OF THE SYMMETRIC EMBEDDED INPUT MATRIX
!
!			= [ INV(M), DEL(INV(M)), DEL(DEL(INV(M))) ]
!
!.....EXTERNAL PROCEDURES
!	EB_LUDCMP	EMBEDDED VERSION OF CHOLESKEY DECOMPOSITION 
!	EB_LUBKSB	EMBEDDED VERSION OF CHOLESKEY BACKSUBSTITUTION
!=================================================
! COPYRIGHT(C) 2003 JAMES D. TURNER, AMDYN SYSTEMS
!=================================================	
USE EB_VARIABLE_HANDLING
IMPLICIT NONE

TYPE(EB),DIMENSION(:,:),INTENT(IN ):: A
TYPE(EB),DIMENSION(:,:),INTENT(OUT):: INV_A 
REAL(DP):: D
INTEGER, DIMENSION(SIZE(A,2))          :: INDX
TYPE(EB),DIMENSION(SIZE(A,2),SIZE(A,2)):: ASAVE, TMP
INTEGER:: I, DIM

DIM = SIZE(A,2)

COPY_A_MATRIX:IF ( SIZE(A,1) /= SIZE(A,2) )THEN
   DO I = 1, DIM
      ASAVE(I,:) = A(I,:)
   END DO !I
ELSE
   ASAVE = A
END IF COPY_A_MATRIX

TMP = EB_EYE(DIM,ONE)

TMP = .constant.0.0D0
DO I = 1, DIM
   TMP(I,I)%E = 1.0D0                                      ! DEFINE IDENTITY MATRIX
END DO

CALL EB_LUDCMP( ASAVE, INDX, D )                      ! LU FACTOR INPUT MATRIX

DO I = 1, DIM                                         ! INVERT MAT ONE COLM AT A TIME
   CALL EB_LUBKSB( ASAVE, INDX, TMP(:,I) )
END DO

COPY_INV_A_MATRIX:IF ( SIZE(A,1) /= SIZE(A,2) )THEN
   INV_A(1:DIM,1:DIM) = TMP(1:DIM,1:DIM)
END IF COPY_INV_A_MATRIX

END SUBROUTINE EB_LU_DECOMP_INV_MATRIX

!==============================================
SUBROUTINE EB_AXeB( DIM, A, b, X_soln )
!  This program inverts embedded linear matrix equations for x of the form:
!
!	Ax = b,  x = inv(A)*b, where inv(A) is obtained by gaussian elimination
!
!    where DIM is the dimension of the linear system
!          A   is the input square matrix consisting of embedded elements
!          b   is the input vector consisting of embedded elements
!          x   is the desired embedded solution vector
!
!   X_soln consists of: [ inv(A)*b, jacobian( inv(A)*b ), hessian( inv(A)*b ) ]
!
!   WHERE TRANSFORMATION TO TRIANGULAR FORM IS GIVEN BY:
!
!   |X X ... X |   |Y|     |X X ... X |   |Y|
!   |X X ... X |   |Y|     |0 X ... X |   |Y|
!   |. . ... . | = |.|  => |0 0 ... X | = |.|
!   |X X ... X |   |Y|     |0 0 ....X |   |Y|  
!
!    NOTE: 
!          No pivoting is used
!          No symmetry is assumed
!	  Set PRT =.TRUE. for detailed solution accuracy printing
!==========================
!copyright (c) 2001 James D. Turner
!==========================  
USE EB_VARIABLE_HANDLING
USE OCEA_IO
IMPLICIT NONE

!.....ARGUMENT LIST VARIABLES

INTEGER,                        INTENT(IN ) ::DIM
TYPE(EB), DIMENSION(dim,dim),   INTENT(IN ) ::A
TYPE(EB), DIMENSION(dim),       INTENT(IN ) ::B
TYPE(EB), DIMENSION(dim),       INTENT(OUT) ::X_SOLN

!.....LOCAL VARIABLES

TYPE(EB), DIMENSION(dim,dim)::A_COPY, A_TRIANG
TYPE(EB), DIMENSION(dim)    ::B_COPY, ERROR, TMP
TYPE(EB)::EB_ZERO, EB_ONE, RECIP_Ajj
INTEGER::I,J,K,N
LOGICAL::PRT,FPASS
character(len=3):: str
SAVE FPASS
SAVE EB_ZERO, EB_ONE
DATA FPASS/.TRUE./

!....Initialize Constants and Matrices

IF (FPASS) THEN
    FPASS    = .FALSE.
    EB_ZERO  = .constant.0.0D0
    EB_ONE   = .constant.0.0D0
    EB_ONE%E = 1.0D0
END IF

N = dim; B_COPY = B; A_TRIANG = A

!	Transform Input embedded Matrix to Triangular Form:  */

DO j= 1,N-1 

   RECIP_Ajj = 1.0D0/A_TRIANG(j,j)

   A_COPY = A_TRIANG

   TMP   (J+1:N) = A_COPY(J+1:N,J)*RECIP_AJJ
   B_COPY(J+1:N) = B_COPY(J+1:N) - TMP(J+1:N)*B_COPY(J)

   DO k=j+1,N
      A_TRIANG(j+1:N,k) = A_COPY(j+1:N,k) - TMP(j+1:N)*A_COPY(j,k) 
   END DO !K

END DO !J

!	Back Substitution for Linear Equation Solution    */

X_SOLN(N) = B_COPY(N)/A_TRIANG(N,N) 

DO J=N-1,1,-1

   RECIP_Ajj = 1.0D0/A_TRIANG(j,j)

   X_SOLN(j) = B_COPY(j)*RECIP_Ajj

   DO k=j+1,N 
      X_SOLN(j) = X_SOLN(j) - ( A_TRIANG(j,k)*RECIP_Ajj )*X_SOLN(k) 
   END DO !K

END DO !J

!...CHECK SOLUTION ACCURACY: AX - B = ERROR

PRT=.FALSE.
CHECK_SOLN_ERROR: IF (PRT) THEN
    ERROR = A*X_SOLN - B
    DO I=1,N
       WRITE(6,*)
       WRITE(6,*)
       write(str,'(i3)')i
       CALL EB_PRT(6,ERROR(I),"AXeB:AX-B=ERROR, I = "//str )
    END DO
END IF CHECK_SOLN_ERROR

END SUBROUTINE EB_AXeB

!===================================
SUBROUTINE EB_MAT_INV_GELIM( DIM, MAT_INPUT, MAT_INV )
!  THIS PROGRAM INVERTS A GENERAL EMBEDDED MATRIX USING GAUSSIAN ELIMINATION.
!  THE ALGORITHM SOLVES FOR ONE COLUMN AT A TIME OF THE INVERSE MATRIX.  THE
!  RESULTING SOLUTION HAS THE FOLLOWING STRUCTURE
!
!		INV(MAT):=[ INV(MAT), JACOBIAN( INV(MAT) ), HESSIAN( INV(MAT) ) ]
!==========================
!copyright (c) 2003 James D. Turner
!==========================  
USE PROBLEM_DATA
USE EB_VARIABLE_HANDLING
IMPLICIT NONE

!...ARGUMENT LIST VARIABLES
INTEGER,                      INTENT(IN ):: DIM
TYPE(EB), DIMENSION(DIM,DIM), INTENT(IN ):: MAT_INPUT
TYPE(EB), DIMENSION(DIM,DIM), INTENT(OUT):: MAT_INV
!...LOCAL VARIABLES
TYPE(EB), DIMENSION(DIM)                 :: IDCOL
INTEGER                                  :: I

! INITIALIZE FIRST COLUMN OF IDENTITY MATRIX COLUMN

IDCOL      = .constant.0.0D0                                 !INITIALIZE ARRAY

IDCOL(1)%E = 1.0D0

! INVERT FOR FIRST COLUMN OF INVERSE MATRIX

CALL EB_AXeB   ( DIM, MAT_INPUT, IDCOL, MAT_INV(:,1) )!GET FIRST COLUMN OF INVERSE MATRIX

! LOOP FOR REMAINING COLUMNES OF INVERSE MATRIX

DO I=2,NV
   IDCOL(I-1)%E = 0.0D0                               !UPDATE RIGHT HAND SIDE FOR AX=B
   IDCOL(I  )%E = 1.0D0
   CALL EB_AXeB( DIM, MAT_INPUT, IDCOL, MAT_INV(:,I) )!GET N-TH  COLUMN OF INVERSE MATRIX
END DO !I

END SUBROUTINE EB_MAT_INV_GELIM

!====================================================
      SUBROUTINE EB_jacobi(a,n,np,d,v,nrot)
      USE EB_VARIABLE_HANDLING
      IMPLICIT NONE

      INTEGER, PARAMETER:: NMAX=500
      INTEGER:: n,np,nrot
      TYPE(EB), DIMENSION(NP,NP), INTENT(INOUT):: a
      TYPE(EB), DIMENSION(NP   ), INTENT(INOUT):: d
      TYPE(EB), DIMENSION(NP,NP), INTENT(INOUT):: v
      INTEGER:: i,ip,iq,j
      TYPE(EB):: c,g,h,s,sm,t,tau,theta,tresh, EB_ONE
      TYPE(EB), DIMENSION(NP):: b, z
      EB_ONE = .constant.0.0D0; EB_ONE%E = ONE
      V = .constant.ZERO
      do 12 ip=1,n
         v(ip,ip)%E=ONE
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=.constant.ZERO
13    continue
      nrot=0
      do 24 i=1,50
        sm=.constant.ZERO
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm%E.eq.ZERO)return
        if(i.lt.4)then
          tresh=.constant.( 0.20D0*sm%E/FLOAT(n**2) )
        else
          tresh=.constant.ZERO
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.0D0*abs(a(ip,iq))
                        if((i.gt.4).and.(abs(d(ip)%E)+g%E.eq.abs(d(ip)%E)) &
                       .and.(abs(d(iq)%E)+g%E.eq.abs(d(iq)%E)))then
              a(ip,iq)=.constant.ZERO
            else if(abs(a(ip,iq)%E).gt.tresh%E)then
              h=d(iq)-d(ip)
              if(abs(h%E)+g%E.eq.abs(h%E))then
                t=a(ip,iq)/h
              else
                theta=0.50D0*h/a(ip,iq)
                t=.constant.ONE/(abs(theta%E)+sqrt(ONE+theta%E**2))
                if(theta%E.lt.ZERO)t=-t
              endif
              c=ONE/sqrt(EB_ONE+t**2)
              s=t*c
              tau=s/(EB_ONE+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=.constant.ZERO
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=.constant.ZERO
23      continue
24    continue
      STOP 'too many iterations in jacobi'
      return
END SUBROUTINE EB_jacobi
!  (C) Copr. 1986-92 Numerical Recipes Software 


!=============================================
FUNCTION EB_EYE(DIM,VALUE)  RESULT(EYE)
    USE EB_VARIABLE_HANDLING
    IMPLICIT NONE

    INTEGER, INTENT(IN):: DIM
    REAL(DP),INTENT(IN):: VALUE
    TYPE(EB), DIMENSION(DIM,DIM):: EYE
    INTEGER:: I
    EYE   = .constant.0.0D0
    DO I=1,DIM
       EYE(I,I)%E = VALUE
    END DO
END FUNCTION EB_EYE


END MODULE OCEA_MATRIX_LIB