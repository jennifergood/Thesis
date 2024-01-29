MODULE OCEA_MATHLIB
USE PROBLEM_DATA
USE EB_VARIABLE_HANDLING
IMPLICIT NONE

!.....DEFINE OPERATOR-OVERLOADED EMBEDDED OBJECT INTRINSIC FUNCTION INTERFACE

   INTERFACE OPERATOR(**)
     MODULE PROCEDURE  EB_RAISED_DP, EB_RAISED_INT
   END INTERFACE

   INTERFACE SIN
        MODULE PROCEDURE EB_SIN
   END INTERFACE

   INTERFACE COS
        MODULE PROCEDURE EB_COS
   END INTERFACE

   INTERFACE TAN
        MODULE PROCEDURE EB_TAN
   END INTERFACE

   INTERFACE ASIN
        MODULE PROCEDURE EB_ASIN
   END INTERFACE

   INTERFACE ACOS
        MODULE PROCEDURE EB_ACOS
   END INTERFACE

   INTERFACE ATAN
        MODULE PROCEDURE EB_ATAN
   END INTERFACE

   INTERFACE ATAN2
        MODULE PROCEDURE EB_EB_ATAN2, EB_DP_ATAN2, DP_EB_ATAN2
   END INTERFACE

   INTERFACE LOG
        MODULE PROCEDURE EB_LN
   END INTERFACE

   INTERFACE EXP
        MODULE PROCEDURE EB_EXP
   END INTERFACE

   INTERFACE SQRT
        MODULE PROCEDURE EB_SQRT
   END INTERFACE

   INTERFACE SINH
        MODULE PROCEDURE EB_SINH
   END INTERFACE

   INTERFACE COSH
        MODULE PROCEDURE EB_COSH
   END INTERFACE

   INTERFACE TANH
        MODULE PROCEDURE EB_TANH
   END INTERFACE

   INTERFACE ABS
        MODULE PROCEDURE EB_ABS, DP_ABS_VEC, DP_ABS
   END INTERFACE
   
   INTERFACE MAG
        MODULE PROCEDURE EB_MAG
   END INTERFACE

   INTERFACE FLOOR
        MODULE PROCEDURE EB_FLOOR
   END INTERFACE

   INTERFACE CEILING
        MODULE PROCEDURE EB_CEILING
   END INTERFACE

!.....DEFINE OVERLOADED EMBEDDED OBJECT OPERATOR FUNCTION INTERFACE

   INTERFACE OPERATOR(.DOT.)   !CONTRACTION OPERATIONS
     MODULE PROCEDURE                                  &
        EB_VEC_DOT_EB_VEC, EB_MAT_DOT_EB_VEC,          &
        DP_VEC_DOT_DP_VEC, DP_MAT_DOT_DP_VEC,          &
        DP_T3_DOT_DP_VEC,  DP_T4_DOT_DP_VEC,           &
        DP_T5_DOT_DP_VEC,  DP_T2_DOT_DP_T3,            &
        DP_T2_DOT_DP_T2,   DP_T3_DOT_DP_T2,            &
        DP_T4_DOT_DP_T2,   DP_T3_DOT_DP_T3,            &
        DP_T2_DOT_DP_T4
   END INTERFACE

   INTERFACE OPERATOR(.T.)
     MODULE PROCEDURE                                  &
        EB_MAT_TRANSPOSE, DP_MAT_TRANSPOSE, EB_VEC_EB_VEC_TRANSPOSE
   END INTERFACE

   INTERFACE OPERATOR(.X.)
     MODULE PROCEDURE                                  &
        EB_VEC_CROSS_EB_VEC, EB_VEC_CROSS_DP_VEC,      &
        DP_VEC_CROSS_EB_VEC, DP_VEC_CROSS_DP_VEC
   END INTERFACE

   INTERFACE OPERATOR(.CONSTANT.)
     MODULE PROCEDURE                                  &
        EB_CONSTANT_EQUALS_DP_SCA
   END INTERFACE


CONTAINS

!...BUILD THE OPERATOR OVERLOADED MATH LIBRARIES

FUNCTION EB_RAISED_DP(A,R)
     TYPE(EB)    EB_RAISED_DP
     TYPE(EB), INTENT(IN):: A
     REAL(DP), INTENT(IN):: R
     REAL(DP), DIMENSION( 0:4):: F

  ! DEFINE FUNCTION AND DERIVATIVES

     F(0) = (A%E)**R

  select case( deriv_order )
     case( 1 )
        F(1) = R*(A%E)**(R-1.0D0)
     case( 2 )
        F(1) =           R*(A%E)**(R-1.0D0)
        F(2) = R*(R-1.0D0)*(A%E)**(R-2.0D0)
     case( 3 )
        F(1) =                     R*(A%E)**(R-1.0D0)
        F(2) =           R*(R-1.0D0)*(A%E)**(R-2.0D0)
        F(3) = R*(R-1.0D0)*(R-2.0D0)*(A%E)**(R-3.0D0)
     case( 4 )
        F(1) =                               R*(A%E)**(R-1.0D0)
        F(2) =                     R*(R-1.0D0)*(A%E)**(R-2.0D0)
        F(3) =           R*(R-1.0D0)*(R-2.0D0)*(A%E)**(R-3.0D0)
        f(4) = R*(R-1.0D0)*(R-2.0D0)*(R-3.0D0)*(A%E)**(R-4.0D0)
  end select

  EB_RAISED_DP = COMPOSITE_FUNCTION( A, F )

END FUNCTION EB_RAISED_DP

FUNCTION EB_RAISED_INT(A,INT)
     TYPE(EB)     EB_RAISED_INT
     TYPE(EB), INTENT(IN):: A
     INTEGER,  INTENT(IN):: INT
     REAL(DP), DIMENSION( 0:4):: F

     ! FUNCTION AND DERIVATIVE DEFINITIONS

  F(0) = (A%E)**INT

  select case( deriv_order )
     case( 1 )
        F(1) = INT*(A%E)**(INT-1)
     case( 2 )
        F(1) =         INT*(A%E)**(INT-1)
        F(2) = INT*(INT-1)*(A%E)**(INT-2)
     case( 3 )
        F(1) =                 INT*(A%E)**(INT-1)
        F(2) =         INT*(INT-1)*(A%E)**(INT-2)
        F(3) = INT*(INT-1)*(INT-2)*(A%E)**(INT-3)
     case( 4 )
        F(1) =                         INT*(A%E)**(INT-1)
        F(2) =                 INT*(INT-1)*(A%E)**(INT-2)
        F(3) =         INT*(INT-1)*(INT-2)*(A%E)**(INT-3)
        f(4) = INT*(INT-1)*(INT-2)*(INT-3)*(A%E)**(INT-4)
  end select

  EB_RAISED_INT = COMPOSITE_FUNCTION( A, F )

END FUNCTION EB_RAISED_INT

FUNCTION EB_SIN( A )
     TYPE(EB)          EB_SIN
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_SIN
     REAL(DP):: S, C
     S = SIN( A%E ); C = COS( A%E )

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_SIN(0) =  S
     CASE(1)
         DERIV_SIN(0) =  S
         DERIV_SIN(1) =  c
     CASE(2)
         DERIV_SIN(0) =  S
         DERIV_SIN(1) =  c
         DERIV_SIN(2) = -S
     CASE(3)
         DERIV_SIN(0) =  S
         DERIV_SIN(1) =  c
         DERIV_SIN(2) = -S
         DERIV_SIN(3) = -C
     CASE(4)
         DERIV_SIN(0) =  S
         DERIV_SIN(1) =  c
         DERIV_SIN(2) = -S
         DERIV_SIN(3) = -C
         DERIV_SIN(4) =  S
     END SELECT 

     EB_SIN = COMPOSITE_FUNCTION( A, DERIV_SIN )

END FUNCTION EB_SIN

FUNCTION EB_COS( A )
     TYPE(EB)          EB_COS
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_COS
     REAL(DP):: S, C
     S = SIN( A%E ); C = COS( A%E )

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_COS(0) =  C
     CASE(1)
         DERIV_COS(0) =  C
         DERIV_COS(1) = -S
     CASE(2)
         DERIV_COS(0) =  C
         DERIV_COS(1) = -S
         DERIV_COS(2) = -C
     CASE(3)
         DERIV_COS(0) =  C
         DERIV_COS(1) = -S
         DERIV_COS(2) = -C
         DERIV_COS(3) =  S
     CASE(4)
         DERIV_COS(0) =  C
         DERIV_COS(1) = -S
         DERIV_COS(2) = -C
         DERIV_COS(3) =  S
         DERIV_COS(4) =  C
     END SELECT 

     EB_COS = COMPOSITE_FUNCTION( A, DERIV_COS )

END FUNCTION EB_COS

FUNCTION EB_TAN( A )
     TYPE(EB)     EB_TAN
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_TAN
     REAL(DP):: S2, T, T2
     T = TAN( A%E ); S2 = ( 1.0D0/COS( A%E ) )**2

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_TAN(0) =  T
     CASE(1)
         DERIV_TAN(0) =  T
         DERIV_TAN(1) =  S2
     CASE(2)
         DERIV_TAN(0) =  T
         DERIV_TAN(1) =  S2
         DERIV_TAN(2) =  2.0D0*S2*T
     CASE(3)
         DERIV_TAN(0) =  T
         DERIV_TAN(1) =  S2
         DERIV_TAN(2) =  2.0D0*S2*T
         DERIV_TAN(3) =  2.0D0*S2*(2.0D0*T*T + S2)
     CASE(4)
         DERIV_TAN(0) =  T
         DERIV_TAN(1) =  S2
         DERIV_TAN(2) =  2.0D0*S2*T
         DERIV_TAN(3) =  2.0D0*S2*(2.0D0*T*T + S2)
         DERIV_TAN(4) =  8.0D0*S2*T*(T*T + 2.0D0*S2)
     END SELECT 

     EB_TAN = COMPOSITE_FUNCTION( A, DERIV_TAN )

END FUNCTION EB_TAN

FUNCTION EB_EXP( A )
     TYPE(EB)     EB_EXP
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_EXP
     REAL(DP):: E
     E = EXP( A%E )

     DERIV_EXP(0:DERIV_ORDER) = E

     EB_EXP = COMPOSITE_FUNCTION( A, DERIV_EXP )

END FUNCTION EB_EXP

FUNCTION EB_LN( A )
     TYPE(EB)     EB_LN
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_LOG
     REAL(DP):: R
     
     IF( A%E  < 0.0D0 ) THEN
         WRITE(6,*)'ERROR:EB_LN A%E < 0 => COMPLEX VARIABLE'
         STOP'OCEA DOES NOT SUPPORT COMPLEX VARIABLES'
     END IF

     R = 1.0D0/A%E

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_LOG(0) =  LOG( A%E )
     CASE(1)
         DERIV_LOG(0) =  LOG( A%E )
         DERIV_LOG(1) =  R
     CASE(2)
         DERIV_LOG(0) =  LOG( A%E )
         DERIV_LOG(1) =  R
         DERIV_LOG(2) = -R*R
     CASE(3)
         DERIV_LOG(0) =  LOG( A%E )
         DERIV_LOG(1) =  R
         DERIV_LOG(2) = -R*R
         DERIV_LOG(3) =  2.0D0*R**3
     CASE(4)
         DERIV_LOG(0) =  LOG( A%E )
         DERIV_LOG(1) =  R
         DERIV_LOG(2) = -R*R
         DERIV_LOG(3) =  2.0D0*R**3
         DERIV_LOG(4) = -6.0D0*R**4
     END SELECT 

     EB_LN = COMPOSITE_FUNCTION( A, DERIV_LOG )

END FUNCTION EB_LN

FUNCTION EB_SQRT( A )
     TYPE(EB)         EB_SQRT
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_SQRT
     REAL(DP):: SQR
     
     IF( A%E  < 0.0D0 ) THEN
         WRITE(6,*)'ERROR:EB_SQRT A%E < 0 => COMPLEX VARIABLE'
         STOP'OCEA DOES NOT SUPPORT COMPLEX VARIABLES'
     END IF
     
     SQR = SQRT( A%E )
     IF ( ABS( SQR ) < 1.0D-13 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR EB_SQRT: A%E = 0 => SINGULAR DERIVATIVES'
         STOP
         END IF

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_SQRT(0) =  SQR
     CASE(1)
         DERIV_SQRT(0) =  SQR
         DERIV_SQRT(1) =  1.0D0/( 2.0D0*SQR   )
     CASE(2)
         DERIV_SQRT(0) =  SQR
         DERIV_SQRT(1) =  1.0D0/( 2.0D0*SQR   )
         DERIV_SQRT(2) = -1.0D0/( 4.0D0*SQR**3)
     CASE(3)
         DERIV_SQRT(0) =  SQR
         DERIV_SQRT(1) =  1.0D0/( 2.0D0*SQR   )
         DERIV_SQRT(2) = -1.0D0/( 4.0D0*SQR**3)
         DERIV_SQRT(3) =  3.0D0/( 8.0D0*SQR**5)
     CASE(4)
         DERIV_SQRT(0) =  SQR
         DERIV_SQRT(1) =  1.0D0/( 2.0D0*SQR   )
         DERIV_SQRT(2) = -1.0D0/( 4.0D0*SQR**3)
         DERIV_SQRT(3) =  3.0D0/( 8.0D0*SQR**5)
         DERIV_SQRT(4) =-15.0D0/(16.0D0*SQR**7)
     END SELECT 

     EB_SQRT = COMPOSITE_FUNCTION( A, DERIV_SQRT )

END FUNCTION EB_SQRT

FUNCTION EB_ACOS( A )
     TYPE(EB)         EB_ACOS
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_ACOS
     REAL(DP):: X, F
     X = A%E; F = 1.0D0/SQRT( 1.0D0 - X*X )

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_ACOS(0) =  ACOS( X )
     CASE(1)
         DERIV_ACOS(0) =  ACOS( X )
         DERIV_ACOS(1) = -F
     CASE(2)
         DERIV_ACOS(0) =  ACOS( X )
         DERIV_ACOS(1) = -F
         DERIV_ACOS(2) = -X*F**3
     CASE(3)
         DERIV_ACOS(0) =  ACOS( X )
         DERIV_ACOS(1) = -F
         DERIV_ACOS(2) = -X*F**3
         DERIV_ACOS(3) = -F**3*( 1.0D0 + 3.0D0*(X*F)**2 )
     CASE(4)
         DERIV_ACOS(0) =  ACOS( X )
         DERIV_ACOS(1) = -F
         DERIV_ACOS(2) = -X*F**3
         DERIV_ACOS(3) = -F**3*( 1.0D0 + 3.0D0*(X*F)**2 )
         DERIV_ACOS(4) = -X*F**5*( 9.0D0 + 15.0D0*(X*F)**2 )
     END SELECT 

     EB_ACOS = COMPOSITE_FUNCTION( A, DERIV_ACOS )

END FUNCTION EB_ACOS

FUNCTION EB_ASIN( A )
     TYPE(EB)         EB_ASIN
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_ASIN
     REAL(DP):: X, F
     X = A%E; F = 1.0D0/SQRT( 1.0D0 - X*X )

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_ASIN(0) =  ASIN( X )
     CASE(1)
         DERIV_ASIN(0) =  ASIN( X )
         DERIV_ASIN(1) =  F
     CASE(2)
         DERIV_ASIN(0) =  ASIN( X )
         DERIV_ASIN(1) =  F
         DERIV_ASIN(2) =  X*F**3
     CASE(3)
         DERIV_ASIN(0) =  ASIN( X )
         DERIV_ASIN(1) =  F
         DERIV_ASIN(2) =  X*F**3
         DERIV_ASIN(3) =  F**3*( 1.0D0 + 3.0D0*(X*F)**2 )
     CASE(4)
         DERIV_ASIN(0) =  ASIN( X )
         DERIV_ASIN(1) =  F
         DERIV_ASIN(2) =  X*F**3
         DERIV_ASIN(3) =  F**3*( 1.0D0 + 3.0D0*(X*F)**2 )
         DERIV_ASIN(4) =  X*F**5*( 9.0D0 + 15.0D0*(X*F)**2 )
     END SELECT 

     EB_ASIN = COMPOSITE_FUNCTION( A, DERIV_ASIN )

END FUNCTION EB_ASIN


FUNCTION EB_ATAN( A )
     TYPE(EB)         EB_ATAN
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_ATAN
     REAL(DP):: X, F
     X = A%E; F = 1.0D0/(1.0D0 + X*X )

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_ATAN(0) =  ATAN( X )
     CASE(1)
         DERIV_ATAN(0) =  ATAN( X )
         DERIV_ATAN(1) =  F
     CASE(2)
         DERIV_ATAN(0) =  ATAN( X )
         DERIV_ATAN(1) =  F
         DERIV_ATAN(2) = -2.0D0*X*F**2
     CASE(3)
         DERIV_ATAN(0) =  ATAN( X )
         DERIV_ATAN(1) =  F
         DERIV_ATAN(2) = -2.0D0*X*F**2
         DERIV_ATAN(3) =  2.0D0*F**2*( 4.0D0*X*X*F - 1.0D0 )
     CASE(4)
         DERIV_ATAN(0) =  ATAN( X )
         DERIV_ATAN(1) =  F
         DERIV_ATAN(2) = -2.0D0*X*F**2
         DERIV_ATAN(3) =  2.0D0*F**2*( 4.0D0*X*X*F - 1.0D0 )
         DERIV_ATAN(4) = 24.0D0*X*F**3*( 1.0D0 - 2.0D0*X*X*F )
     END SELECT 

     EB_ATAN = COMPOSITE_FUNCTION( A, DERIV_ATAN )

END FUNCTION EB_ATAN


FUNCTION EB_EB_ATAN2( N, D )
     TYPE(EB)         EB_EB_ATAN2
     TYPE(EB), INTENT(IN)    :: N
     TYPE(EB), INTENT(IN)    :: D
     REAL(DP), DIMENSION(0:4):: DERIV_ATAN2
     REAL(DP)                :: X, F
! 
     X = N%E/D%E; F = 1.0D0/( 1.0D0 + X*X )

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_ATAN2(0) =  ATAN2( N%E, D%E )
     CASE(1)
         DERIV_ATAN2(0) =  ATAN2( N%E, D%E )
         DERIV_ATAN2(1) =  F
     CASE(2)
         DERIV_ATAN2(0) =  ATAN2( N%E, D%E )
         DERIV_ATAN2(1) =  F
         DERIV_ATAN2(2) = -2.0D0*X*F**2
     CASE(3)
         DERIV_ATAN2(0) =  ATAN2( N%E, D%E )
         DERIV_ATAN2(1) =  F
         DERIV_ATAN2(2) = -2.0D0*X*F**2
         DERIV_ATAN2(3) =  2.0D0*F**2*( 4.0D0*X*X*F - 1.0D0 )
     CASE(4)
         DERIV_ATAN2(0) =  ATAN2( N%E, D%E )
         DERIV_ATAN2(1) =  F
         DERIV_ATAN2(2) = -2.0D0*X*F**2
         DERIV_ATAN2(3) =  2.0D0*F**2*( 4.0D0*X*X*F - 1.0D0 )
         DERIV_ATAN2(4) = 24.0D0*X*F**3*( 1.0D0 - 2.0D0*X*X*F )
     END SELECT 

     EB_EB_ATAN2 = COMPOSITE_FUNCTION( N/D, DERIV_ATAN2 )

END FUNCTION EB_EB_ATAN2


FUNCTION DP_EB_ATAN2( N, D )
     TYPE(EB)        DP_EB_ATAN2
     REAL(DP), INTENT(IN):: N
     TYPE(EB), INTENT(IN):: D
     REAL(DP), DIMENSION(0:4):: DERIV_ATAN2
     REAL(DP):: X, F

     X = N/D%E; F = 1.0D0/( 1.0D0 + X*X )

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_ATAN2(0) =  ATAN2( N, D%E )
     CASE(1)
         DERIV_ATAN2(0) =  ATAN2( N, D%E )
         DERIV_ATAN2(1) =  F
     CASE(2)
         DERIV_ATAN2(0) =  ATAN2( N, D%E )
         DERIV_ATAN2(1) =  F
         DERIV_ATAN2(2) = -2.0D0*X*F**2
     CASE(3)
         DERIV_ATAN2(0) =  ATAN2( N, D%E )
         DERIV_ATAN2(1) =  F
         DERIV_ATAN2(2) = -2.0D0*X*F**2
         DERIV_ATAN2(3) =  2.0D0*F**2*( 4.0D0*X*X*F - 1.0D0 )
     CASE(4)
         DERIV_ATAN2(0) =  ATAN2( N, D%E )
         DERIV_ATAN2(1) =  F
         DERIV_ATAN2(2) = -2.0D0*X*F**2
         DERIV_ATAN2(3) =  2.0D0*F**2*( 4.0D0*X*X*F - 1.0D0 )
         DERIV_ATAN2(4) = 24.0D0*X*F**3*( 1.0D0 - 2.0D0*X*X*F )
     END SELECT 

     DP_EB_ATAN2 = COMPOSITE_FUNCTION( N/D, DERIV_ATAN2 )

END FUNCTION DP_EB_ATAN2


FUNCTION EB_DP_ATAN2( N, D )
     TYPE(EB)        EB_DP_ATAN2
     TYPE(EB), INTENT(IN):: N
     REAL(DP), INTENT(IN):: D
     REAL(DP), DIMENSION(0:4):: DERIV_ATAN2
     REAL(DP):: X, F
!
     X = N%E/D; F = 1.0D0/( 1.0D0 + X*X )

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_ATAN2(0) =  ATAN2( N%E, D )
     CASE(1)
         DERIV_ATAN2(0) =  ATAN2( N%E, D )
         DERIV_ATAN2(1) =  F
     CASE(2)
         DERIV_ATAN2(0) =  ATAN2( N%E, D )
         DERIV_ATAN2(1) =  F
         DERIV_ATAN2(2) = -2.0D0*N%E*X*F**2
     CASE(3)
         DERIV_ATAN2(0) =  ATAN2( N%E, D )
         DERIV_ATAN2(1) =  F
         DERIV_ATAN2(2) = -2.0D0*X*F**2
         DERIV_ATAN2(3) =  2.0D0*F**2*( 4.0D0*X*X*F - 1.0D0 )
     CASE(4)
         DERIV_ATAN2(0) =  ATAN2( N%E, D )
         DERIV_ATAN2(1) =  F
         DERIV_ATAN2(2) = -2.0D0*X*F**2
         DERIV_ATAN2(3) =  2.0D0*F**2*( 4.0D0*X*X*F - 1.0D0 )
         DERIV_ATAN2(4) = 24.0D0*X*F**3*( 1.0D0 - 2.0D0*X*X*F )
     END SELECT 

     EB_DP_ATAN2 = COMPOSITE_FUNCTION( N/D, DERIV_ATAN2 )

END FUNCTION EB_DP_ATAN2

FUNCTION EB_SINH( A )
     TYPE(EB)          EB_SINH
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_SINH
     REAL(DP):: S, C
     S = SINH( A%E ); C = COSH( A%E )

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_SINH(0) =  S
     CASE(1)
         DERIV_SINH(0) =  S
         DERIV_SINH(1) =  c
     CASE(2)
         DERIV_SINH(0) =  S
         DERIV_SINH(1) =  c
         DERIV_SINH(2) =  S
     CASE(3)
         DERIV_SINH(0) =  S
         DERIV_SINH(1) =  c
         DERIV_SINH(2) =  S
         DERIV_SINH(3) =  C
     CASE(4)
         DERIV_SINH(0) =  S
         DERIV_SINH(1) =  c
         DERIV_SINH(2) =  S
         DERIV_SINH(3) =  C
         DERIV_SINH(4) =  S
     END SELECT 

     EB_SINH = COMPOSITE_FUNCTION( A, DERIV_SINH )

END FUNCTION EB_SINH

FUNCTION EB_COSH( A )
     TYPE(EB)          EB_COSH
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_COSH
     REAL(DP):: S, C
     S = SINH( A%E ); C = COSH( A%E )

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_COSH(0) =  C
     CASE(1)
         DERIV_COSH(0) =  C
         DERIV_COSH(1) =  S
     CASE(2)
         DERIV_COSH(0) =  C
         DERIV_COSH(1) =  S
         DERIV_COSH(2) =  C
     CASE(3)
         DERIV_COSH(0) =  C
         DERIV_COSH(1) =  S
         DERIV_COSH(2) =  C
         DERIV_COSH(3) =  S
     CASE(4)
         DERIV_COSH(0) =  C
         DERIV_COSH(1) =  S
         DERIV_COSH(2) =  C
         DERIV_COSH(3) =  S
         DERIV_COSH(4) =  C
     END SELECT 

     EB_COSH = COMPOSITE_FUNCTION( A, DERIV_COSH )

END FUNCTION EB_COSH

FUNCTION EB_TANH( A )
     TYPE(EB)     EB_TANH
     TYPE(EB), INTENT(IN):: A
     REAL(DP), DIMENSION(0:4):: DERIV_TANH
     REAL(DP):: S2, T, T2
     T = TANH( A%E ); S2 = ( 1.0D0/COSH( A%E ) )**2

     SELECT CASE (DERIV_ORDER)
     CASE(0)
         DERIV_TANH(0) =  T
     CASE(1)
         DERIV_TANH(0) =  T
         DERIV_TANH(1) =  S2
     CASE(2)
         DERIV_TANH(0) =  T
         DERIV_TANH(1) =  S2
         DERIV_TANH(2) = -2.0D0*S2*T
     CASE(3)
         DERIV_TANH(0) =  T
         DERIV_TANH(1) =  S2
         DERIV_TANH(2) = -2.0D0*S2*T
         DERIV_TANH(3) =  2.0D0*S2*(2.0D0*T*T - S2)
     CASE(4)
         DERIV_TANH(0) =  T
         DERIV_TANH(1) =  S2
         DERIV_TANH(2) = -2.0D0*S2*T
         DERIV_TANH(3) =  2.0D0*S2*(2.0D0*T*T - S2)
         DERIV_TANH(4) =  8.0D0*S2*T*(2.0D0*S2 - T*T)
     END SELECT 

     EB_TANH = COMPOSITE_FUNCTION( A, DERIV_TANH )

END FUNCTION EB_TANH

FUNCTION EB_ABS( A )
     TYPE(EB)          EB_ABS
     TYPE(EB), INTENT(IN):: A
     REAL(DP)            :: F, DF
     ! FUNCTION AND DERIVATIVE DEFINITIONS
       F = ABS( A%E )                       ! FUNCTION
      DF = SIGN( 1.0D0, A%E )               ! FIRST DERIVATIVE => ABS(X)/X
     ! EMBEDDED FUNCTION DEFINITIONS
     EB_ABS%E  = F
     EB_ABS%V  = DF*A%V
     EB_ABS%T2 = DF*A%T2
     EB_ABS%T3 = DF*A%T3
     EB_ABS%T4 = DF*A%T4
END FUNCTION EB_ABS

FUNCTION DP_ABS_VEC( A )
     REAL(DP)                        :: DP_ABS_VEC
     TYPE(EB),DIMENSION(:),INTENT(IN):: A
     
     DP_ABS_VEC = SQRT( A%E .DOT. A%E )
     
END FUNCTION DP_ABS_VEC

FUNCTION DP_ABS( A )
     REAL(DP)                        :: DP_ABS
     REAL(DP),DIMENSION(:),INTENT(IN):: A
     
     DP_ABS = SQRT( A(1)*A(1) + A(2)*A(2) + A(3)*A(3) )
     
END FUNCTION DP_ABS

FUNCTION EB_MAG( A )            RESULT( MAGNITUDE )
     TYPE(EB),DIMENSION(:),INTENT(IN):: A
     REAL(DP)                        :: MAGNITUDE
     
     MAGNITUDE = SQRT( A%E .DOT. A%E ) 
!                VECTOR NORM

END FUNCTION EB_MAG

FUNCTION EB_FLOOR( A )
     TYPE(EB)        EB_FLOOR
     TYPE(EB), INTENT(IN):: A
     REAL(DP)            :: F
     ! FUNCTION AND DERIVATIVE DEFINITIONS
       F = FLOOR( A%E )
     ! EMBEDDED FUNCTION DEFINITIONS
     EB_FLOOR%E  = F
     EB_FLOOR%V  = 0.0D0
     EB_FLOOR%T2 = 0.0D0
     EB_FLOOR%T3 = 0.0D0
     EB_FLOOR%T4 = 0.0D0
END FUNCTION EB_FLOOR

FUNCTION EB_CEILING( A )
     TYPE(EB)      EB_CEILING
     TYPE(EB), INTENT(IN):: A
     REAL(DP)            :: F
     ! FUNCTION AND DERIVATIVE DEFINITIONS
       F = FLOOR( A%E )
     ! EMBEDDED FUNCTION DEFINITIONS
     EB_CEILING%E  = F
     EB_CEILING%V  = 0.0D0
     EB_CEILING%T2 = 0.0D0
     EB_CEILING%T3 = 0.0D0
     EB_CEILING%T4 = 0.0D0
END FUNCTION EB_CEILING


!...UTILITIES FOR OVERLOADING THE ".DOT.= DOT PRODUCT" ASSIGNMENT OPERATOR

FUNCTION  DP_VEC_DOT_DP_VEC(A,B)  RESULT(VdotV)
     REAL(DP)::                          VdotV
     REAL(DP), DIMENSION(:), INTENT(IN):: A,B
     ! EMBEDDED FUNCTION DEFINITIONS
     VdotV = DOT_PRODUCT( A, B )
END FUNCTION DP_VEC_DOT_DP_VEC

FUNCTION DP_MAT_DOT_DP_VEC(A,B)      RESULT(MdotV)
     REAL(DP), DIMENSION(:,:), INTENT(IN)::   A
     REAL(DP), DIMENSION(:  ), INTENT(IN)::   B
     REAL(DP), DIMENSION(SIZE(B))        :: MdotV
     INTEGER:: I,J
     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-VECTOR MULTIPLICATION

     MdotV = MATMUL( A, B )

END FUNCTION DP_MAT_DOT_DP_VEC

FUNCTION DP_T3_DOT_DP_VEC(A,B)           RESULT(T3dotV)
     ! DP MATRIX ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:),   INTENT(IN)::    A
     REAL(DP), DIMENSION(:    ),   INTENT(IN)::    B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)):: T3dotV
     INTEGER:: I,J,K

     IF( SIZE(A,2) /= SIZE(B) .OR. SIZE(A,3) /= SIZE(B) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T3_DOT_DP_VEC INCONSISTENT MATRIX-VECTOR DIMS'
        WRITE(6,*)'CURRENT MATRIX DIMS =',SIZE(A,1),SIZE(A,2),SIZE(A,3)
        WRITE(6,*)'CURRENT VECTOR DIM  =',SIZE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-VECTOR MULTIPLICATION
     
     DO I = 1, SIZE(B)
        DO J = 1, SIZE(B)
           T3dotV(I,J) = DOT_PRODUCT( A(I,J,1:SIZE(B)), B(1:SIZE(B)) )
        END DO !J
     END DO !I

END FUNCTION DP_T3_DOT_DP_VEC

FUNCTION DP_T4_DOT_DP_VEC(A,B)     RESULT( T4dotV )
     ! DP MATRIX ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:,:),   INTENT(IN)::    A
     REAL(DP), DIMENSION(:    ),     INTENT(IN)::    B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(B),SIZE(B)):: T4dotV
 
     INTEGER::I,J,K,L
     
     DO I = 1, SIZE(B)
        DO J = 1, SIZE(B)
           DO K = 1, SIZE(B)
              T4dotV(I,J,K) = DOT_PRODUCT( A(I,J,K, 1:SIZE(B)), B(1:SIZE(B)) ) 
           END DO !K
        END DO ! J
     END DO ! I

END FUNCTION DP_T4_DOT_DP_VEC

FUNCTION DP_T5_DOT_DP_VEC(A,B)     RESULT( T5dotV )
     ! DP MATRIX ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:,:,:), INTENT(IN)::    A
     REAL(DP), DIMENSION(:    ),     INTENT(IN)::    B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(B),SIZE(B),SIZE(B)):: T5dotV
 
     INTEGER::I,J,K,L
     
     DO I = 1, SIZE(B)
        DO J = 1, SIZE(B)
           DO K = 1, SIZE(B)
              DO L = 1, SIZE(B)
              T5dotV(I,J,K,L) = DOT_PRODUCT( A(I,J,K,L,1:SIZE(B)), B(1:SIZE(B)) )
              END DO !L
           END DO !K
        END DO ! J
     END DO ! I

END FUNCTION DP_T5_DOT_DP_VEC

FUNCTION EB_VEC_DOT_EB_VEC(A,B)    RESULT(VdotV)
     TYPE(EB)                             VdotV
     TYPE(EB), DIMENSION(:), INTENT(IN):: A,B
     INTEGER:: I

     IF( SIZE(A) /= SIZE(B) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:EB_VEC_DOT_EB_VEC DIM A/=B'
        WRITE(6,*)'CURRENT VECTOR A DIM  =',SIZE(A)
        WRITE(6,*)'CURRENT VECTOR B DIM  =',SIZE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     VdotV = A(1)*B(1)
     DO I = 2, SIZE(B)
        VdotV = VdotV + A(I)*B(I)
     END DO

END FUNCTION EB_VEC_DOT_EB_VEC

FUNCTION EB_MAT_DOT_EB_VEC(A,B)      RESULT(MdotV)
     ! EB VECTOR ARRAY FUNCTION
     TYPE(EB), DIMENSION(:,:), INTENT(IN)::   A
     TYPE(EB), DIMENSION(:  ), INTENT(IN)::   B
     TYPE(EB), DIMENSION(SIZE(A))        :: MdotV
     INTEGER:: I,J

     IF( SIZE(A,2) /= SIZE(B) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:EB_MAT_DOT_EB_VEC DIM A(*,2) /= DIM B'
        WRITE(6,*)'CURRENT MATRIX DIM =',SIZE(A,1),SIZE(A,2)
        WRITE(6,*)'CURRENT VECTOR DIM =',SIZE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-VECTOR MULTIPLICATION

     DO I=1,SIZE(A,1)
        MdotV(I) = A(I,1)*B(1)
        DO J=2,SIZE(A,2)
           MdotV(I) = MdotV(I) + A(I,j)*B(j)
        END DO
     END DO

END FUNCTION EB_MAT_DOT_EB_VEC

FUNCTION DP_T2_DOT_DP_T2(A,B)      RESULT(T2dotT2)
     ! EB VECTOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:     ), INTENT(IN)::   A
     REAL(DP), DIMENSION(:,:     ), INTENT(IN)::   B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(B,2)) :: T2dotT2
     INTEGER:: n1,n2,s

     IF( SIZE(A,2) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T2_DOT_DP_T2 DIM A(*,2) /= DIM B'
        WRITE(6,*)'CURRENT MATRIX DIM =',SIZE(A,1),SIZE(A,2)
        WRITE(6,*)'CURRENT VECTOR DIM =',SIZE(B,1)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-TENSOR MULTIPLICATION
     
     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(B,2)
           T2dotT2(n1,n2) = dot_product( A(n1,:), B(:,n2) )
        END DO
     END DO

END FUNCTION DP_T2_DOT_DP_T2

FUNCTION DP_T2_DOT_DP_T3(A,B)      RESULT(T2dotT3)
     ! EB VECTOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:    ), INTENT(IN)::   A
     REAL(DP), DIMENSION(:,:,:  ), INTENT(IN)::   B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(B,2),SIZE(B,3)):: T2dotT3
     INTEGER:: n1,n2,n3,s

     IF( SIZE(A,2) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T2_DOT_DP_T3 DIM A(*,2) /= DIM B'
        WRITE(6,*)'CURRENT MATRIX DIM =',SIZE(A,1),SIZE(A,2)
        WRITE(6,*)'CURRENT VECTOR DIM =',SIZE(B,1),SIZE(B,2),SIZE(B,3)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-TENSOR MULTIPLICATION
     
     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(B,2)
           DO n3=1,SIZE(B,3)
              T2dotT3(n1,n2,n3) = dot_product( A(n1,:), B(:,n2,n3) )
           END DO
        END DO
     END DO

END FUNCTION DP_T2_DOT_DP_T3

FUNCTION DP_T3_DOT_DP_T2(A,B)      RESULT(T3dotT2)
     ! EB VECTOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:    ), INTENT(IN)::   A
     REAL(DP), DIMENSION(:,:      ), INTENT(IN)::   B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(B,2)):: T3dotT2
     INTEGER:: n1,n2,n3,s

     IF( SIZE(A,3) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T3_DOT_DP_T2 DIM A(*,2) /= DIM B'
        WRITE(6,*)'CURRENT MATRIX DIM =',SIZE(A,1),SIZE(A,2),SIZE(A,3)
        WRITE(6,*)'CURRENT VECTOR DIM =',SIZE(B,1),SIZE(B,2)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-TENSOR MULTIPLICATION
     
     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(A,2)
           DO n3=1,SIZE(B,2)
              T3dotT2(n1,n2,n3) = dot_product( A(n1,n2,:), B(:,n3) )
           END DO
        END DO
     END DO

END FUNCTION DP_T3_DOT_DP_T2

FUNCTION DP_T4_DOT_DP_T2(A,B)      RESULT(T4dotT2)
     ! EB VECTOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:,:), INTENT(IN)::   A
     REAL(DP), DIMENSION(:,:    ), INTENT(IN)::   B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(A,3),SIZE(B,2)):: T4dotT2
     INTEGER:: n1,n2,n3,n4,s

     IF( SIZE(A,4) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T4_DOT_DP_T2 DIM A(*,2) /= DIM B'
        WRITE(6,*)'CURRENT MATRIX DIM =',SIZE(A,1),SIZE(A,2),SIZE(A,3),SIZE(A,4)
        WRITE(6,*)'CURRENT VECTOR DIM =',SIZE(B,1),SIZE(B,2)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-TENSOR MULTIPLICATION A(I,J,K,:)*B(:,M)
     
     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(A,2)
           DO n3=1,SIZE(A,3)
              DO n4=1,SIZE(B,2)
                 T4dotT2(n1,n2,n3,n4) = dot_product( A(n1,n2,n3,:), B(:,n4) )
              END DO
           END DO
        END DO
     END DO

END FUNCTION DP_T4_DOT_DP_T2

FUNCTION DP_T2_DOT_DP_T4(A,B)      RESULT(T2dotT4)
     ! EB VECTOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:      ), INTENT(IN)::   A
     REAL(DP), DIMENSION(:,:,:,:  ), INTENT(IN)::   B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(B,2),SIZE(B,3),SIZE(B,4)):: T2dotT4
     INTEGER:: n1,n2,n3,n4,s

     IF( SIZE(A,2) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T2_DOT_DP_T4 DIM A(*,2) /= DIM B'
        WRITE(6,*)'CURRENT MATRIX DIM =',SIZE(A,1),SIZE(A,2)
        WRITE(6,*)'CURRENT VECTOR DIM =',SIZE(B,1),SIZE(B,2),SIZE(B,3),SIZE(B,4)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-TENSOR MULTIPLICATION A(I,:)*B(:,K,L,M)
     
     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(B,2)
           DO n3=1,SIZE(B,3)
              DO n4=1,SIZE(B,4)
                 T2dotT4(n1,n2,n3,n4) = dot_product( A(n1,:), B(:,n2,n3,n4) )
              END DO
           END DO
        END DO
     END DO

END FUNCTION DP_T2_DOT_DP_T4

FUNCTION DP_T3_DOT_DP_T3(A,B)      RESULT(T3dotT3)
     ! EB VECTOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:    ), INTENT(IN)::   A
     REAL(DP), DIMENSION(:,:,:    ), INTENT(IN)::   B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(B,2),SIZE(B,3)):: T3dotT3
     INTEGER:: n1,n2,n3,n4,s

     IF( SIZE(A,3) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T3_DOT_DP_T3 DIM A(*,2) /= DIM B'
        WRITE(6,*)'CURRENT MATRIX DIM =',SIZE(A,1),SIZE(A,2),SIZE(A,3)
        WRITE(6,*)'CURRENT VECTOR DIM =',SIZE(B,1),SIZE(B,2),SIZE(B,3)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-TENSOR MULTIPLICATION A(I,J,:)*B(:,K,L)
         
     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(A,2)
           DO n3=1,SIZE(B,2)
              DO n4=1,SIZE(B,3)
                 T3dotT3(n1,n2,n3,n4) = dot_product( a(n1,n2,:), b(:,n3,n4) )
              END DO
           END DO
        END DO
     END DO

END FUNCTION DP_T3_DOT_DP_T3

FUNCTION DP_MAT_TRANSPOSE(A)               RESULT(TMAT)
     REAL(DP), DIMENSION(:,:), INTENT(IN):: A
     REAL(DP), DIMENSION( SIZE(A,2), SIZE(A,1) )::TMAT

     ! EMBEDDED FUNCTION DEFINITIONS

     TMAT = TRANSPOSE( A )

END FUNCTION DP_MAT_TRANSPOSE

FUNCTION EB_MAT_TRANSPOSE(A)               RESULT(TMAT)
     ! EB MATRIX ARRAY FUNCTION
     TYPE(EB), DIMENSION(:,:), INTENT(IN):: A
     TYPE(EB), DIMENSION( SIZE(A,2), SIZE(A,1) )::TMAT
     INTEGER:: J,N,M
     N = SIZE(A,1); M = SIZE(A,2)

     ! EMBEDDED FUNCTION DEFINITIONS

     DO J = 1, M
        TMAT(J,1:N) = A(1:N,J)
     END DO

END FUNCTION EB_MAT_TRANSPOSE

FUNCTION EB_VEC_EB_VEC_TRANSPOSE( A, B )  RESULT( OUTTER )
    TYPE(EB), DIMENSION(:), INTENT(IN):: A,B
    TYPE(EB), DIMENSION(SIZE(A),SIZE(B)):: OUTTER
    INTEGER:: I,J

       DO J = 1, SIZE(B)
          OUTTER(:,J) = A(:)*B(J)
       END DO

END FUNCTION EB_VEC_EB_VEC_TRANSPOSE

!...UTILITIES FOR OVERLOADING THE ".X.= CROSS PRODUCT" ASSIGNMENT OPERATOR

FUNCTION EB_VEC_CROSS_EB_VEC(A,B) RESULT(AxB)
    ! EB VECTOR ARRAY FUNCTION
    TYPE(EB), DIMENSION(3), INTENT(IN):: A,B
    TYPE(EB), DIMENSION(3)            :: AxB

    ! EMBEDDED VECTOR CROSS PRODUCT OPERATOR

    AxB(1) = A(2)*B(3) - A(3)*B(2)
    AxB(2) = A(3)*B(1) - A(1)*B(3)
    AxB(3) = A(1)*B(2) - A(2)*B(1)

END FUNCTION EB_VEC_CROSS_EB_VEC

FUNCTION EB_VEC_CROSS_DP_VEC(A,B) RESULT(AxB)
    ! EB VECTOR ARRAY FUNCTION
    TYPE(EB), DIMENSION(3), INTENT(IN):: A
    REAL(DP), DIMENSION(3), INTENT(IN):: B
    TYPE(EB), DIMENSION(3)            :: AxB

    ! EMBEDDED VECTOR CROSS PRODUCT OPERATOR

    AxB(1) = A(2)*B(3) - A(3)*B(2)
    AxB(2) = A(3)*B(1) - A(1)*B(3)
    AxB(3) = A(1)*B(2) - A(2)*B(1)

END FUNCTION EB_VEC_CROSS_DP_VEC

FUNCTION DP_VEC_CROSS_EB_VEC(A,B) RESULT(AxB)
    ! EB VECTOR ARRAY FUNCTION
    REAL(DP), DIMENSION(3), INTENT(IN):: A
    TYPE(EB), DIMENSION(3), INTENT(IN):: B
    TYPE(EB), DIMENSION(3)            :: AxB

    ! EMBEDDED VECTOR CROSS PRODUCT OPERATOR

    AxB(1) = A(2)*B(3) - A(3)*B(2)
    AxB(2) = A(3)*B(1) - A(1)*B(3)
    AxB(3) = A(1)*B(2) - A(2)*B(1)

END FUNCTION DP_VEC_CROSS_EB_VEC

FUNCTION DP_VEC_CROSS_DP_VEC(A,B) RESULT(AxB)
    ! EB VECTOR ARRAY FUNCTION
    REAL(DP), DIMENSION(3), INTENT(IN):: A,B
    REAL(DP), DIMENSION(3)            :: AxB

    ! EMBEDDED VECTOR CROSS PRODUCT OPERATOR

    AxB(1) = A(2)*B(3) - A(3)*B(2)
    AxB(2) = A(3)*B(1) - A(1)*B(3)
    AxB(3) = A(1)*B(2) - A(2)*B(1)

END FUNCTION DP_VEC_CROSS_DP_VEC

!....DEFINE EMBEDDED OBJECT CONSTANT

FUNCTION EB_CONSTANT_EQUALS_DP_SCA(A) RESULT(EB_CONSTANT)
    REAL(DP),INTENT(IN):: A
    TYPE(EB)           :: EB_CONSTANT
    
    SELECT CASE (DERIV_ORDER)
      CASE(0)
        EB_CONSTANT%E = A
      CASE(1)
        EB_CONSTANT%E = A
        EB_CONSTANT%V = 0.0D0
      CASE(2)
        EB_CONSTANT%E  = A
        EB_CONSTANT%V  = 0.0D0
        EB_CONSTANT%T2 = 0.0D0
      CASE(3)
        EB_CONSTANT%E  = A
        EB_CONSTANT%V  = 0.0D0
        EB_CONSTANT%T2 = 0.0D0
        EB_CONSTANT%T3 = 0.0D0
      CASE(4)
        EB_CONSTANT%E  = A
        EB_CONSTANT%V  = 0.0D0
        EB_CONSTANT%T2 = 0.0D0
        EB_CONSTANT%T3 = 0.0D0
        EB_CONSTANT%T4 = 0.0D0
    END SELECT

END FUNCTION EB_CONSTANT_EQUALS_DP_SCA


!     ==================
!.....COMPOSITE FUNCTION.....
!     ================== 

function COMPOSITE_FUNCTION( Y, F )   RESULT( Z )! COMPOSITE FUNCTION    
!  INPUT
!   Y   OCEA VARIABLE
!	F	FUNCTION AND DERIV VALUES
!=================================
!COPYRIGHT (C) 2003 JAMES D. TURNER
!==================================  
    type(eb),intent(in)               :: y
    real(dp),dimension(NV)            :: y1
    real(dp),dimension(NV,NV)         :: y2
    real(dp),dimension(NV,NV,NV)      :: y3
    real(dp),dimension(NV,NV,NV,NV)   :: y4
    real(dp),DIMENSION(0:4),intent(in):: f
    TYPE(EB):: Z
    integer:: j,k,r

!.....UNPACK OCEA VARIALBE

    SELECT CASE( DERIV_ORDER )
      CASE(1)
         Y1(:) = Y%v%VPART(:)
      CASE(2)
         Y1(:)   = Y%v%VPART(:)
         Y2(:,:) = Y%t2%T2PART(:,:)
      CASE(3)
         Y1(:)     = Y%v%vpart(:)
         Y2(:,:)   = Y%t2%T2PART(:,:)
         Y3(:,:,:) = Y%t3%T3PART(:,:,:)
      CASE(4)
         Y1(:)       = Y%v%VPART(:)
         Y2(:,:)     = Y%t2%T2PART(:,:)
         Y3(:,:,:)   = Y%t3%T3PART(:,:,:)
         Y4(:,:,:,:) = Y%t4%T4PART(:,:,:,:)
    END SELECT


!.....PACK OCEA COMPOSITE FUNCTION SUB-OBJECTS

    SELECT CASE( DERIV_ORDER )
      CASE(0)

         Z%E = F(0)

      CASE(1)

         Z%E = F(0)
         Z%V%vpart(:) = F(1)*Y1(:)

      CASE(2)

         Z%E = F(0)
         Z%V%vpart(:)        = F(1)*Y1(:)
         DO J = 1, NV
            Z%T2%T2PART(:,J) = F(2)*Y1(:)*Y1(J) + F(1)*Y2(:,J)
         END DO !J

      CASE(3)

         Z%E                      = F(0)
         Z%V%vpart(:)             = F(1)*Y1(:)
         DO J = 1, NV
            Z%T2%T2PART(:,J)      = F(2)*Y1(:)*Y1(J) + F(1)*Y2(:,J)
            DO K = 1, NV
               Z%T3%T3PART(:,J,K) =                                           &
                 F(3)*  Y1(:)*Y1(J)*Y1(K)                               +     &
                 F(2)*( Y2(:,K)*Y1(J) + Y1(:)*Y2(J,K) + Y2(:,J)*Y1(K) ) +     &
                 F(1)*Y3(:,J,K)
            END DO !K
         END DO !J

      CASE(4)

         Z%E                      = F(0)
         Z%V%vpart(:)             = F(1)*Y1(:)
         DO J = 1, NV
            Z%T2%T2PART(:,J)      = F(2)*Y1(:)*Y1(J) + F(1)*Y2(:,J)
            DO K = 1, NV
               Z%T3%T3PART(:,J,K) =                                           &
                 F(3)*  Y1(:)*Y1(J)*Y1(K)                               +     &
                 F(2)*( Y2(:,K)*Y1(J) + Y1(:)*Y2(J,K) + Y2(:,J)*Y1(K) ) +     &
                 F(1)*  Y3(:,J,K)
               DO R = 1, NV
                  Z%T4%T4PART(:,J,K,R) =                                      &
                    F(4)*  Y1(:)*Y1(J)*Y1(K)*Y1(R)                      +     &
                    F(3)*( Y2(:,R)*Y1(J)*Y1(K) + Y1(:)*Y2(J,R)*Y1(K)    +     &
                           Y1(:)*Y1(J)*Y2(K,R) + Y2(:,K)*Y1(J)*Y1(R)    +     &
                           Y1(:)*Y2(J,K)*Y1(R) + Y2(:,J)*Y1(K)*Y1(R) )  +     &
                    F(2)*( Y3(:,K,R)*Y1(J)     + Y2(:,R)*Y2(J,K)        +     &
                           Y3(:,J,R)*Y1(K)     + Y2(:,K)*Y2(J,R)        +     &
                           Y1(:)*Y3(J,K,R)     + Y2(:,J)*Y2(K,R)        +     &
                           Y3(:,J,K)*Y1(R)                           )  +     &
                    F(1)*  Y4(:,J,K,R)
               END DO !R
            END DO !K
         END DO !J

      END SELECT

END FUNCTION COMPOSITE_FUNCTION

      
FUNCTION BINOMIAL( M, N ) RESULT( BINOM )
!==========================================================
!
!__PURPOSE:	
!       This routine computes the standard binomial coefficient
!	for integer variables M and N.  The Basic equation
!	being computed is given by:
!
!		BINOM = M!/N!/(M-N)!
!
!.....CAN HANDLE 400 DOF FOR 4TH ORDER EXPANSION (VERIFIED WITH MACSYMA)
!.....CAN HANDLE 100 DOF FOR 5TH ORDER EXPANSION (VERIFIED WITH MACSYMA)
!.....CAN HANDLE 100 DOF FOR 6TH ORDER EXPANSION (VERIFIED WITH MACSYMA)
!
!.....ISSUE: THE MAXIMUM SIZE INTEGER *4 VARIABLE IS: 2147483647 ON A 32-BIT COMPUTER
!
!==========================================================
!__COPYRIGHT (c) 1997 James D. Turner
!  CONVERTED TO F90 2006 JAMES D. TURNER
!==========================================================
	implicit none
	integer,              INTENT(IN):: M, N
	INTEGER                         :: BINOM
	INTEGER,ALLOCATABLE,DIMENSION(:):: NUM, DOM
	INTEGER                         :: N_NUM, N_DOM
	INTEGER                         :: BIG_DOWN
    INTEGER                         :: I, J, SMALL
    
      SMALL = MIN( M, N )
    
      IF(      SMALL == 0 )THEN
          BINOM = 1; RETURN
      ELSE IF( SMALL == 1 ) THEN
          BINOM = M; RETURN
      ELSE IF( SMALL > 0 .AND. M == N ) THEN
          BINOM = 1; RETURN
      END IF
!.....QUICK RETURN CHECK    
!    
	  BIG_DOWN = MAX( N, M - N )
	    
	  N_NUM    = M - BIG_DOWN
	  N_DOM    = MIN( N, M - N ) - 1
!.....COMPUTE FACTORS TO MINIMIZE CALCULATIONS
	    
	  ALLOCATE( NUM( N_NUM ) )
	  ALLOCATE( DOM( N_DOM ) )
	    
	  FORALL( I = 1:N_NUM ) NUM(I) = BIG_DOWN + I
	    
	  FORALL( I = 1:N_DOM ) DOM(I) = I + 1
!.....ALLOCATE AND FILL ARRAYS WITH REMAINING FACTORIAL FACTORS	    

	    DO I = N_DOM, 1, -1
	       DO J = 1, N_NUM
	          IF( .NOT. DOM(I) == 1 .AND. MOD( NUM(J), DOM(I) ) == 0 ) THEN
	              NUM(J) = NUM(J)/DOM(I)
	              DOM(I) = 1
	          END IF
	       END DO ! J
	       IF( SUM(DOM) == N_DOM ) EXIT
	    END DO !I
!.....ELIMINATE COMMON FACTORS BEFORE FORMING PRODUCT
      
      IF( SUM(DOM) /= N_DOM ) THEN
         WRITE(6,*)'BINOMIAL:ERROR SUM(DOM) /= N_DOM REDUCTION FAILED, N,M=',N,M
      END IF
	    
	  BINOM = PRODUCT( NUM )
!.....FORM BINOMIAL FROM REDUCED FACTORS

      IF( BINOM < 0 ) THEN
         WRITE(6,*)'BINOMIAL:ERROR BINOM < 0 OUT OF RANGE CALCULATION N, M =',N,M
         STOP
      END IF

	  DEALLOCATE( NUM )
	  DEALLOCATE( DOM )

	end FUNCTION BINOMIAL

END MODULE OCEA_MATHLIB