MODULE EB_VARIABLE_HANDLING
USE PROBLEM_DATA
IMPLICIT NONE

!NAMING CONVENTION FOR OVERLOADED FUNCTIONS AND SUBROUTINES:
!
!  ( DATA TYPE  ) ( ARRAY TYPE ) OPERATOR ( DATA TYPE ) ( ARRAY TYPE )
!
!    DATA TYPE:  
!        DP => DOUBLE PRECISION
!		 EB => OCEA DERIVED DATA TYPE( DOUBLE PRECISION )
!        IT => INTEGER
!    ARRAY TYPE: SCA=> SCALAR
!		 VEC=> Nx1 VECTOR ARRAY
!		 MAT=> NxN MATRIX ARRAY
!
!    OCEA_ARRAY:= { T1, T2, T3, T4, T5, T6 } Nx1 EMBEDDED DATA OBJECT
!
!	 T1 => Nx1         DP ARRAY (FETCH OCEA DATA T1=OCEA_ARRAY)
!	 T2 => NxN         DP ARRAY (FETCH OCEA DATA T2=OCEA_ARRAY)
!	 T3 => NxNxN       DP ARRAY (FETCH OCEA DATA T3=OCEA ARRAY)
!	 T4 => NxNxNxN     DP ARRAY (FETCH OCEA DATA T4=OCEA ARRAY)
!    T5 => NxNxNxNxN   DP ARRAY (FETCH OCEA DATA T5=OCEA ARRAY)
!    T6 => NxNxNxNxNxN DP ARRAY (FETCH OCEA DATA T6=OCEA ARRAY)
!
!    NOTE:
!          F:={ F, DF, DDF, DDDF, DDDDF }
!
!   F=SCALAR:={ F, T1, T2,  T3,   T4    }  => T3 = F  OVERLOADING FOR DATA FETCHING
!
!   F=VECTOR:={ T1,T2, T3,  T4,   T5    }  => T4 = F  OVERLOADING FOR DATA FETCHING
!
!   F=MATRIX:={ T2,T3, T4,  T5,   T6    }  => T6 = F  OVERLOADING FOR DATA FETCHING

!.....DEFINE USER-DEFINED DATA TYPES

   TYPE PV
       REAL(DP) VPART(NV)
   END TYPE PV

   TYPE PT2
       REAL(DP) T2PART(NV,NV)
   END TYPE PT2

   TYPE PT3
       REAL(DP) T3PART(NV,NV,NV)
   END TYPE PT3
   
   TYPE PT4
       REAL(DP) T4PART(NV,NV,NV,NV)
   END TYPE PT4

   TYPE PT5
       REAL(DP) T5PART(NV,NV,NV,NV,NV)
   END TYPE PT5

   TYPE EB
       REAL(DP ) E           ! 1x1       scalar
       TYPE(PV ) V           ! nx1       vector
       TYPE(PT2) T2          ! nxn       matrix
       TYPE(PT3) T3          ! nxnxn     tensor
       TYPE(PT4) T4          ! nxnxnxn   tensor
       TYPE(PT5) T5          ! nxnxnxnxn tensor
   END TYPE EB

!.....DEFINE INTRINSIC OPERATOR INTERFACES FOR OVERLOADING

   INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE  &
        PV_EQ_DP,       & ! pv(nx1)         = dp(1x1) vector eqn
       PT2_EQ_DP,       & !pt2(nxn)         = dp(1x1) matrix eqn
       PT3_EQ_DP,       & !pt3(nxnxn)       = dp(1x1) tensor eqn
       PT4_EQ_DP,       & !pt4(nxnxnxn)     = dp(1x1) tensor eqn
       PT5_EQ_DP,       & !pt5(nxnxnxn)     = dp(1x1) tensor eqn
       
       DP_EQ_PV_SCA,    & ! DP(1X1)         =PV (1X1) SUB-OBJECT TYPE FETCHING
       DP_EQ_PT2_SCA,   & ! DP(1X1)         =PT2(1X1) SUB-OBJECT TYPE FETCHING
       DP_EQ_PT3_SCA,   & ! DP(1X1)         =PT3(1X1) SUB-OBJECT TYPE FETCHING
       DP_EQ_PT4_SCA,   & ! DP(1X1)         =PT4(1X1) SUB-OBJECT TYPE FETCHING
       DP_EQ_PT5_SCA,   & ! DP(1X1)         =PT5(1X1) SUB-OBJECT TYPE FETCHING
                        
       E_EQ_E,          & !  e(1x1)         =  e(1x1) scalar eqn
       E_EQ_DP,         & !  e(1x1)         = dp(1x1) scalar eqn
        
       DP_T1_EQ_EB_SCA, & ! dp(nx1)         =  e(1x1) fetch vector part of scalar eqn
       DP_T2_EQ_EB_SCA, & ! dp(nxn)         =  e(1x1) fetch matrix part of scalar eqn
       DP_T3_EQ_EB_SCA, & ! dp(nxnxn)       =  e(1x1) fetch tensor part of scalar eqn
       DP_T4_EQ_EB_SCA, & ! dp(nxnxnxn)     =  e(1x1) fetch tensor part of scalar eqn
       DP_T5_EQ_EB_SCA, & ! dp(nxnxnxnxn)   =  e(1x1) fetch tensor part of scalar eqn
       
       DP_T1_EQ_EB_VEC, & ! dp(nx1)         =  e(nx1) fetch vector part of vector eqn
       DP_T2_EQ_EB_VEC, & ! dp(nxn)         =  e(nx1) fetch matrix part of vector eqn
       DP_T3_EQ_EB_VEC, & ! dp(nxnxn)       =  e(nx1) fetch tensor part of vector eqn
       DP_T4_EQ_EB_VEC, & ! dp(nxnxnxn)     =  e(nx1) fetch tensor part of vector eqn
       DP_T5_EQ_EB_VEC, & ! dp(nxnxnxnxn)   =  e(nx1) fetch tensor part of vector eqn
                                        
       DP_T2_EQ_EB_MAT, & ! dp(nxn)         =  e(nxn) fetch matrix part of matrix eqn
       DP_T3_EQ_EB_MAT, & ! dp(nxnxn)       =  e(nxn) fetch tensor part of matrix eqn
       DP_T4_EQ_EB_MAT, & ! dp(nxnxnxn)     =  e(nxn) fetch tensor part of matrix eqn
       DP_T5_EQ_EB_MAT, & ! dp(nxnxnxnxn)   =  e(nxn) fetch tensor part of matrix eqn
       DP_T6_EQ_EB_MAT    ! dp(nxnxnxnxnxn) =  e(nxn) fetch tensor part of matrix eqn

   END INTERFACE

   INTERFACE OPERATOR(+)
      MODULE PROCEDURE  &
       PV_ADD_PV, PT2_ADD_Pt2, PT3_ADD_Pt3, PT4_ADD_Pt4, &!SUB-OBJ OVERLOADING                      
        E_ADD_E, E_ADD_DP, DP_ADD_E                       !FUNC OVERLOADING
   END INTERFACE

   INTERFACE OPERATOR(-)
      MODULE PROCEDURE  &
       PV_SUB_PV, PT2_SUB_Pt2, PT3_SUB_Pt3, PT4_SUB_Pt4,  &!SUB-OBJ OVERLOADING
        E_SUB_E,    MINUS_E,   E_MINUS_DP,  DP_MINUS_E,   &!FUNC OVERLOADING 
        MINUS_PV,   MINUS_PT2,   MINUS_PT3,    MINUS_PT4   !SUB-OBJECT OVERLOADING                       
   END INTERFACE

   INTERFACE OPERATOR(*)
      MODULE PROCEDURE  &
       PV_MULT_DP,      & !pv(1x1)  *  dp(1x1) scalar eqn
      PT2_MULT_DP,      & !pt2(1x1) *  dp(1x1) scalar eqn
      PT3_MULT_DP,      & !pt3(1x1) *  dp(1x1) scalar eqn
      PT4_MULT_DP,      & !pt4(1x1) *  dp(1x1) scalar eqn
      
       DP_MULT_PV,      & !dp(1x1)  *  pv(1x1) scalar eqn
       DP_MULT_PT2,     & !dp(1x1)  * pt2(1x1) scalar eqn
       DP_MULT_PT3,     & !dp(1x1)  * pt3(1x1) scalar eqn
       DP_MULT_PT4,     & !dp(1x1)  * pt4(1x1) scalar eqn   
                            
        E_MULT_E,       & ! e(1x1)  *   e(1x1) scalar eqn
        E_MULT_DP,      & ! e(1x1)  *  dp(1x1) scalar eqn
       DP_MULT_E,       & !dp(1x1)  *   e(1x1) scalar eqn
        
   EB_VEC_MULT_DP_SCA,  & ! e(nx1)  *  dp(1x1) vector eqn
   DP_SCA_MULT_EB_VEC,  & !dp(1x1)  *   e(nx1) vector eqn
   EB_VEC_MULT_EB_SCA,  & ! e(nx1)  *   e(1x1) vector eqn
   EB_SCA_MULT_EB_VEC,  & ! e(1x1)  *   e(nx1) vector eqn
   
   EB_MAT_MULT_DP_SCA,  & ! e(nxn)  *  dp(1x1) matrix eqn
   DP_SCA_MULT_EB_MAT,  & !dp(1x1)  *   e(nxn) matrix eqn
   EB_MAT_MULT_EB_SCA,  & ! e(nxn)  *   e(1x1) matrix eqn
   EB_SCA_MULT_EB_MAT,  & ! e(1x1)  *   e(nxn) matrix eqn
 
   DP_MAT_MULT_EB_VEC,  & !dp(nxn)  *   e(nx1) vector eqn
   EB_MAT_MULT_DP_VEC,  & ! e(nxn)  *  dp(nx1) vector eqn
   EB_MAT_MULT_EB_VEC,  & ! e(nxn)  *   e(nx1) vector eqn
   
   EB_MAT_MULT_EB_MAT,  & ! e(nxn)  *   e(nxn) matrix eqn
   EB_MAT_MULT_DP_MAT,  & ! e(nxn)  *  dp(nxn) matrix eqn
   DP_MAT_MULT_EB_MAT     !dp(nxn)  *   e(nxn) matrix eqn
   
   END INTERFACE

   INTERFACE OPERATOR(/)
      MODULE PROCEDURE &
        E_DIV_E, E_DIV_DP, DP_DIV_E                        !FUNC OVERLOADING
   END INTERFACE
   


CONTAINS

!.....DEFINE SUB-OBJECT ASSIGNMENT RULES FOR (=)

elemental subroutine pv_eq_dp(A,B)
   type( pv ), intent(out)::a
   real( dp ), intent(in )::b
   a%vpart(:) = b
end subroutine pv_eq_dp

elemental subroutine pT2_eq_dp(A,B)
   type( pt2 ), intent(out)::a
   real( dp  ), intent(in )::b
   a%t2part(:,:) = b
end subroutine pT2_eq_dp

elemental subroutine pt3_eq_dp(A,B)
   type( pt3 ), intent(out)::a
   real( dp  ), intent(in )::b
   a%t3part(:,:,:) = b
end subroutine pt3_eq_dp

elemental subroutine pt4_eq_dp(A,B)
   type( pt4 ), intent(out)::a
   real( dp  ), intent(in )::b
   a%t4part(:,:,:,:) = b
end subroutine pt4_eq_dp

elemental subroutine pt5_eq_dp(A,B)
   type( pt5 ), intent(out)::a
   real( dp  ), intent(in )::b
   a%t5part(:,:,:,:,:) = b
end subroutine pt5_eq_dp
!.....DEFINE OCEA OPERATORS FOR (=)

SUBROUTINE DP_EQ_PV_SCA(A,B)
   REAL( DP  ), INTENT(OUT):: A
   TYPE( PV  ), INTENT(IN ):: B
   A = B%VPART(1)
END SUBROUTINE DP_EQ_PV_SCA

SUBROUTINE DP_EQ_PT2_SCA(A,B)
   REAL( DP  ), INTENT(OUT):: A
   TYPE( PT2 ), INTENT(IN ):: B
   A = B%T2PART(1,1)
END SUBROUTINE DP_EQ_PT2_SCA

SUBROUTINE DP_EQ_PT3_SCA(A,B)
   REAL( DP  ), INTENT(OUT):: A
   TYPE( PT3 ), INTENT(IN ):: B
   A = B%T3PART(1,1,1)
END SUBROUTINE DP_EQ_PT3_SCA

SUBROUTINE DP_EQ_PT4_SCA(A,B)
   REAL( DP  ), INTENT(OUT):: A
   TYPE( PT4 ), INTENT(IN ):: B
   A = B%T4PART(1,1,1,1)
END SUBROUTINE DP_EQ_PT4_SCA

SUBROUTINE DP_EQ_PT5_SCA(A,B)
   REAL( DP  ), INTENT(OUT):: A
   TYPE( PT5 ), INTENT(IN ):: B
   A = B%T5PART(1,1,1,1,1)
END SUBROUTINE DP_EQ_PT5_SCA

elemental subroutine e_eq_e(A,B)
   USE PROBLEM_DATA
   type( eb ),intent(out):: A
   type( eb ),intent(in ):: B

   select case( deriv_order )
      case(0)
        A%E = B%E
      case(1)
        A%E = B%E; A%V = B%V
      CASE(2)
        A%E = B%E; A%V = B%V; A%T2 = B%T2
      CASE(3)
        A%E = B%E; A%V = B%V; A%T2 = B%T2; A%T3 = B%T3
      CASE(4)
        A%E = B%E; A%V = B%V; A%T2 = B%T2; A%T3 = B%T3; A%T4 = B%T4
      END SELECT 

end subroutine e_eq_e

ELEMENTAL SUBROUTINE E_EQ_DP(A,B)
   TYPE( EB ),INTENT(OUT):: A
   REAL( DP ),INTENT(IN ):: B
   
   select case( deriv_order )
      case(0)
        A%E = B
      case(1)
        A%E = B; A%V = B
      CASE(2)
        A%E = B; A%V = B; A%T2 = B
      CASE(3)
        A%E = B; A%V = B; A%T2 = B; A%T3 = B
      CASE(4)
        A%E = B; A%V = B; A%T2 = B; A%T3 = B; A%T4 = B
      CASE(5)
        A%E = B; A%V = B; A%T2 = B; A%T3 = B; A%T4 = B; A%T5 = B
   END SELECT 
   
END SUBROUTINE E_EQ_DP

!.....UTILITIES FOR EXTRACTING PARTS OF AN EMBEDDED SCALAR IN TERMS OF REAL ARRAYS

SUBROUTINE DP_T1_EQ_EB_SCA(A,B) ! FETCHS JACOBIAN OF A SCALAR
     REAL(DP), DIMENSION(:), INTENT(OUT):: A
     TYPE(EB),               INTENT(IN ):: B

     ! EMBEDDED FUNCTION DEFINITIONS

     A(:) = B%V%VPART(:)
     
END SUBROUTINE DP_T1_EQ_EB_SCA

SUBROUTINE DP_T2_EQ_EB_SCA(A,B) ! FETCHS HESSIAN OF A SCALAR
     REAL(DP), DIMENSION(:,:), INTENT(OUT):: A
     TYPE(EB),                 INTENT(IN ):: B

     ! CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 2 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T2_EQ_EB_SCA--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 2 FOR SCALAR INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 2'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     A(:,:) = B%T2%T2PART(:,:)
     
END SUBROUTINE DP_T2_EQ_EB_SCA

SUBROUTINE DP_T3_EQ_EB_SCA(A,B) ! FETCHS JACOBIAN( HESSIAN ) OF A SCALAR
     REAL(DP), DIMENSION(:,:,:), INTENT(OUT):: A
     TYPE(EB),                   INTENT(IN ):: B

     ! CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 3 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T3_EQ_EB_SCA--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 3 FOR SCALAR INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 3'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     A(:,:,:) = B%T3%T3PART(:,:,:)

END SUBROUTINE DP_T3_EQ_EB_SCA

SUBROUTINE DP_T4_EQ_EB_SCA(A,B) ! FETCHS HESSIAN( HESSIAN ) OF A SCALAR
     REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT):: A
     TYPE(EB),                     INTENT(IN ):: B

     ! CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 4 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T4_EQ_EB_SCA--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 4 FOR SCALAR INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 4'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     A(:,:,:,:) = B%T4%T4PART(:,:,:,:)

END SUBROUTINE DP_T4_EQ_EB_SCA

SUBROUTINE DP_T5_EQ_EB_SCA(A,B) ! FETCHS HESSIAN( HESSIAN ) OF A SCALAR
     REAL(DP), DIMENSION(:,:,:,:,:), INTENT(OUT):: A
     TYPE(EB),                       INTENT(IN ):: B

     ! CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 5 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T5_EQ_EB_SCA--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 5 FOR SCALAR INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 5'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     A(:,:,:,:,:) = B%T5%T5PART(:,:,:,:,:)

END SUBROUTINE DP_T5_EQ_EB_SCA

!.....UTILITIES FOR EXTRACTING PARTS OF AN EMBEDDED VECTOR IN TERMS OF REAL ARRAYS

SUBROUTINE DP_T1_EQ_EB_VEC(A,B) !FETCHS 0TH ORDER DERIVATIVE OF A VECTOR
     ! DP VECTOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:), INTENT(OUT):: A
     TYPE(EB), DIMENSION(:), INTENT(IN ):: B

     !BEGIN DIMENSIONAL ERROR CHECKING
     
     IF( SIZE(A) /= SIZE(B) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T1_EQ_EB_VEC INPUT VECTORS NOT THE SAME LENGTH'
        WRITE(6,*)'ARRAYS A AND B MUST HAVE SAME LEADING DIM'
        WRITE(6,*)'DIM(A), DIM(B) =',SIZE(A),SIZE(B)
        STOP'INCONSISTENT VECTOR LENGHS ON INPUT'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     A = B%E

END SUBROUTINE DP_T1_EQ_EB_VEC

SUBROUTINE DP_T2_EQ_EB_VEC(A,B) !FETCHS 1ST ORDER DERIVATIVE OF A VECTOR
     ! DP MATRIX ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:), INTENT(OUT):: A
     TYPE(EB), DIMENSION(:),   INTENT(IN ):: B
     INTEGER:: I,J

!.....CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 1 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T2_EQ_EB_VEC--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 1 FOR VECTOR INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 1'
     END IF

!.....BEGIN DIMENSIONAL ERROR CHECKING

     DO I = 1, 2
        IF( SIZE(A,I) /= NV ) THEN
            WRITE(6,*)
            WRITE(6,*)'ERROR:DP_T2_EQ_EB_VEC OUTPUT DIMENSION ERROR DETECTED'
            DO J = 1, 2
               WRITE(6,*)'ITH DIMENSION =',SIZE(A,J)
            END DO
            STOP'INCONSISTENT DIMENSIONS FOR OUTPUT'
        END IF
     END DO !I

!.....EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, NV
        A(:,I) = B%V%VPART(I)
     END DO !I

END SUBROUTINE DP_T2_EQ_EB_VEC

SUBROUTINE DP_T3_EQ_EB_VEC(A,B) !FETCHS 2ND ORDER DERIVATIVE OF A VECTOR
     ! DP 3ED-ORDER TENSOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:), INTENT(OUT):: A
     TYPE(EB), DIMENSION(:),     INTENT(IN ):: B
     INTEGER:: I,J

!.....CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 2 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T3_EQ_EB_VEC--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 2 FOR VECTOR INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 2'
     END IF

!.....BEGIN DIMENSIONAL ERROR CHECKING

     DO I = 1, 3
        IF( SIZE(A,I) /= NV ) THEN
            WRITE(6,*)
            WRITE(6,*)'ERROR:DP_T3_EQ_EB_VEC OUTPUT DIMENSION ERROR DETECTED'
            DO J = 1, 3
               WRITE(6,*)'ITH DIMENSION =',SIZE(A,J)
            END DO
            STOP'INCONSISTENT DIMENSIONS FOR OUTPUT'
        END IF
     END DO !I

!.....EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, NV
        DO J = 1, NV
           A(:,I,J) = B%T2%T2PART(I,J)
        END DO !J
     END DO !I

END SUBROUTINE DP_T3_EQ_EB_VEC

SUBROUTINE DP_T4_EQ_EB_VEC(A,B) !FETCHS 3ED ORDER DERIVATIVE OF A VECTOR
     ! DP 4ED-ORDER TENSOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT):: A
     TYPE(EB), DIMENSION(:),       INTENT(IN ):: B
     INTEGER:: I,J,K

!.....CHECK DERIVATIVE ORDER

     !IF( DERIV_ORDER < 3 .OR. DERIV_ORDER > MAX_DERIV_ORDER ) THEN
     !    WRITE(6,*)
     !    IF( DERIV_ORDER < 3 ) THEN
     !        WRITE(6,*)'ERROR:DP_T4_EQ_EB_VEC--UNDEFINED OPERATION'
     !        WRITE(6,*)'DERIV_ORDER MUST BE >= 3 FOR VECTOR INPUTS'
     !       WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
     !        END IF
     !    IF( DERIV_ORDER > MAX_DERIV_ORDER ) THEN
     !        WRITE(6,*)'ERROR:DP_T4_EQ_EB_VEC--UNDEFINED OPERATION'
     !        WRITE(6,*)'DERIV_ORDER > MAX_DERIVATIVE_ORDER'
     !       WRITE(6,*)'DERIV_ORDER MUST BE < MAX_DERIVATIVE_ORDER FOR VECTOR INPUTS'
     !       WRITE(6,*)'MAX_DERIVATIVE_ORDER =',MAX_DERIV_ORDER
     !        END IF
     !    STOP 'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 3'
     !END IF

!.....BEGIN DIMENSIONAL ERROR CHECKING

     DO I = 1, 4
        IF( SIZE(A,I) /= NV ) THEN
            WRITE(6,*)
            WRITE(6,*)'ERROR:DP_T4_EQ_EB_VEC OUTPUT DIMENSION ERROR DETECTED'
            DO J = 1, 4
               WRITE(6,*)'ITH DIMENSION =',SIZE(A,J)
            END DO
            STOP'INCONSISTENT DIMENSIONS FOR OUTPUT'
        END IF
     END DO !I

!.....EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, NV
        DO J = 1, NV
           DO K = 1, NV
              A(:,I,J,K) = B%T3%T3PART(I,J,K)
           END DO !K
        END DO !J
     END DO !I

END SUBROUTINE DP_T4_EQ_EB_VEC

SUBROUTINE DP_T5_EQ_EB_VEC(A,B) !FETCHS 4TH ORDER DERIVATIVE OF A VECTOR
     ! DP 5ED-ORDER TENSOR ARRAY FUNCTION
     USE PROBLEM_DATA
     REAL(DP), DIMENSION(:,:,:,:,:), INTENT(OUT):: A
     TYPE(EB), DIMENSION(:),         INTENT(IN ):: B
     INTEGER:: I,J,K,L

!.....CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 4 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T5_EQ_EB_VEC--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 4 FOR VECTOR INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 4'
     END IF

!.....BEGIN DIMENSIONAL ERROR CHECKING

     DO I = 1, 5
        IF( SIZE(A,I) /= NV ) THEN
            WRITE(6,*)
            WRITE(6,*)'ERROR:DP_T5_EQ_EB_VEC OUTPUT DIMENSION ERROR DETECTED'
            DO J = 1, 5
               WRITE(6,*)'ITH DIMENSION =',SIZE(A,J)
            END DO
            STOP'INCONSISTENT DIMENSIONS FOR OUTPUT'
        END IF
     END DO !I

!.....EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, NV
        DO J = 1, NV
           DO K = 1, NV
              DO L = 1, NV
                 A(:,I,J,K,L) = B%T4%T4PART(I,J,K,L)
              END DO !L
           END DO !K
        END DO !J
     END DO !I

END SUBROUTINE DP_T5_EQ_EB_VEC

!.....UTILITIES FOR EXTRACTING PARTS OF AN EMBEDDED MATRIX IN TERMS OF REAL ARRAYS

SUBROUTINE DP_T2_EQ_EB_MAT(A,B) !FETCHS 0TH ORDER DERIVATIVE OF A MATRIX
     ! DP MATRIX ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:), INTENT(OUT):: A
     TYPE(EB), DIMENSION(:,:), INTENT(IN ):: B
     INTEGER:: I,J

!.....BEGIN DIMENSIONAL ERROR CHECKING

     DO I = 1, 2
        IF( SIZE(A,I) /= SIZE(B,I) ) THEN
            WRITE(6,*)
            WRITE(6,*)'ERROR:DP_T2_EQ_EB_MAT OUTPUT DIMENSION ERROR DETECTED'
            DO J = 1, 2
               WRITE(6,*)'ITH DIMENSION =',SIZE(A,J),SIZE(B,J)
            END DO
            STOP'INCONSISTENT DIMENSIONS FOR OUTPUT'
        END IF
     END DO !I

     ! EMBEDDED FUNCTION DEFINITIONS

     A = B%E

END SUBROUTINE DP_T2_EQ_EB_MAT

SUBROUTINE DP_T3_EQ_EB_MAT(A,B) !FETCHS 1ST ORDER DERIVATIVE OF A MATRIX
   USE PROBLEM_DATA
     ! DP 3ED-ORDER TENSOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:), INTENT(OUT):: A
     TYPE(EB), DIMENSION(:,:),   INTENT(IN ):: B
     INTEGER:: I,J

!.....CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 1 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T3_EQ_EB_MAT--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 1 FOR MATRIX INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 1'
     END IF

!.....BEGIN DIMENSIONAL ERROR CHECKING

     DO I = 1, 2
        IF( SIZE(A,I) /= SIZE(B,I) .OR. SIZE(A,3) /= NV ) THEN
            WRITE(6,*)
            WRITE(6,*)'ERROR:DP_T3_EQ_EB_MAT OUTPUT DIMENSION ERROR DETECTED'
            DO J = 1, 2
               WRITE(6,*)'ITH DIMENSION =',SIZE(A,J),SIZE(B,J)
            END DO
            WRITE(6,*)'3ED DIMENSION =',SIZE(A,3)
            STOP'INCONSISTENT DIMENSIONS FOR OUTPUT'
        END IF
     END DO !I

     ! EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, NV
        A(:,:,I) = B%V%VPART(I)
     END DO !I

END SUBROUTINE DP_T3_EQ_EB_MAT

SUBROUTINE DP_T4_EQ_EB_MAT(A,B) !FETCHS 2ND ORDER DERIVATIVE OF A MATRIX
   USE PROBLEM_DATA
     ! DP 4TH-ORDER TENSOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT):: A
     TYPE(EB), DIMENSION(:,:),     INTENT(IN ):: B
     INTEGER:: I,J

!.....CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 2 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T4_EQ_EB_MAT--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 2 FOR MATRIX INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 2'
     END IF

!.....BEGIN DIMENSIONAL ERROR CHECKING

     DO I = 1, 2
        IF( SIZE(A,I) /= SIZE(B,I) .OR. SIZE(A,I+2) /= NV ) THEN
            WRITE(6,*)
            WRITE(6,*)'ERROR:DP_T4_EQ_EB_MAT OUTPUT DIMENSION ERROR DETECTED'
            DO J = 1, 2
               WRITE(6,*)'ITH DIMENSION =',SIZE(A,J),SIZE(B,J)
            END DO
            WRITE(6,*)'3ED DIMENSION =',SIZE(A,3)
            WRITE(6,*)'4TH DIMENSION =',SIZE(A,4)
            STOP'INCONSISTENT DIMENSIONS FOR OUTPUT'
        END IF
     END DO !I

!.....EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, NV
        DO J = 1, NV
           A(:,:,I,J) = B%T2%T2PART(I,J)
        END DO !J
     END DO !I

END SUBROUTINE DP_T4_EQ_EB_MAT

SUBROUTINE DP_T5_EQ_EB_MAT(A,B) !FETCHS 3ED ORDER DERIVATIVE OF A MATRIX
    USE PROBLEM_DATA
    ! DP 5TH-ORDER TENSOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:,:,:), INTENT(OUT):: A
     TYPE(EB), DIMENSION(:,:),       INTENT(IN ):: B
     INTEGER:: I,J,K

!.....CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 3 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T5_EQ_EB_MAT--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 3 FOR MATRIX INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 3'
     END IF

!.....BEGIN DIMENSIONAL ERROR CHECKING

     DO I = 1, 2
        IF( SIZE(A,I  ) /= SIZE(B,I) .OR. &
            SIZE(A,I+2) /= NV        .OR. &
            SIZE(A,5  ) /= NV               ) THEN
            WRITE(6,*)
            WRITE(6,*)'ERROR:DP_T5_EQ_EB_MAT OUTPUT DIMENSION ERROR DETECTED'
            DO J = 1, 2
               WRITE(6,*)'ITH DIMENSION =',SIZE(A,J),SIZE(B,J)
            END DO
            WRITE(6,*)'3ED DIMENSION =',SIZE(A,3)
            WRITE(6,*)'4TH DIMENSION =',SIZE(A,4)
            WRITE(6,*)'5TH DIMENSION =',SIZE(A,5)
            STOP'INCONSISTENT DIMENSIONS FOR OUTPUT'
        END IF
     END DO !I

!.....EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, NV
        DO J = 1, NV
           DO K = 1, NV
              A(:,:,I,J,K) = B%T3%T3PART(I,J,K)
           END DO !K
        END DO !J
     END DO !I

END SUBROUTINE DP_T5_EQ_EB_MAT

SUBROUTINE DP_T6_EQ_EB_MAT(A,B) !FETCHS 4TH ORDER DERIVATIVE OF A MATRIX
   USE PROBLEM_DATA
     ! DP 6TH-ORDER TENSOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:,:,:,:,:), INTENT(OUT):: A
     TYPE(EB), DIMENSION(:,:),         INTENT(IN ):: B
     INTEGER:: I,J,K,L

!.....CHECK DERIVATIVE ORDER

     IF( DERIV_ORDER < 4 ) THEN
         WRITE(6,*)
         WRITE(6,*)'ERROR:DP_T6_EQ_EB_MAT--UNDEFINED OPERATION'
         WRITE(6,*)'DERIV_ORDER MUST BE >= 4 FOR MATRIX INPUTS'
         WRITE(6,*)'DERIV_ORDER =',DERIV_ORDER
         STOP'DATA FETCHING OPERATION UNDEFINED: DERIV_ORDER < 4'
     END IF

!.....BEGIN DIMENSIONAL ERROR CHECKING

     DO I = 1, 2
        IF( SIZE(A,I) /= SIZE(B,I) .OR. SIZE(A,I+2) /= NV ) THEN
            WRITE(6,*)
            WRITE(6,*)'ERROR:DP_T4_EQ_EB_MAT OUTPUT DIMENSION ERROR DETECTED'
            DO J = 1, 2
               WRITE(6,*)'ITH DIMENSION =',SIZE(A,J),SIZE(B,J)
            END DO
            WRITE(6,*)'3ED DIMENSION =',SIZE(A,3)
            WRITE(6,*)'4TH DIMENSION =',SIZE(A,4)
            STOP'INCONSISTENT DIMENSIONS FOR OUTPUT'
        END IF
     END DO !I

!.....EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, NV
        DO J = 1, NV
           DO K = 1, NV
              DO L = 1, NV
                 A(:,:,I,J,K,L) = B%T4%T4PART(I,J,K,L)
              END DO !L
           END DO !K
        END DO !J
     END DO !I

END SUBROUTINE DP_T6_EQ_EB_MAT

!.....DEFINE SUB-OPERATOR RULES FOR (+)

elemental function pv_add_pv(A,B)
   type( pv ), intent(in )::a,b
   type( pv )             ::pv_add_pv
   integer::i
   do i = 1, nv
      pv_add_pv%vpart(I) = a%vpart(I) + b%vpart(I)
   end do
end function pv_add_pv

elemental function pt2_add_pt2(A,B)
   type( pt2 ), intent(in )::a,b
   type( pt2 )             ::pt2_add_pt2
   integer::i
   do i = 1, nv
      pt2_add_pt2%t2part(:,I) = a%t2part(:,I) + b%t2part(:,I)
   end do
end function pt2_add_pt2

elemental function pt3_add_pt3(A,B)
   type( pt3 ), intent(in )::a,b
   type( pt3 )             ::pt3_add_pt3
   integer::i
   do i = 1, nv
      pt3_add_pt3%t3part(:,:,I) = a%t3part(:,:,I) + b%t3part(:,:,I)
   end do
end function pt3_add_pt3

elemental function pt4_add_pt4(A,B)
   type( pt4 ), intent(in )::a,b
   type( pt4 )             ::pt4_add_pt4
   integer::i
   do i = 1, nv
      pt4_add_pt4%t4part(:,:,:,I) = a%t4part(:,:,:,I) + b%t4part(:,:,:,I)
   end do
end function pt4_add_pt4

!.....DEFINE OCEA OPERATOR RULES FOR (+)

elemental function e_add_e(A,B)  
   USE PROBLEM_DATA
   type( eb ),intent(in )::A,B
   type( eb )            ::e_add_e

   select case( deriv_order )
      case(0)
        e_add_e%E  = A%E  + B%E
      case(1)
        e_add_e%E  = A%E  + B%E
        e_add_e%V  = A%V  + B%V
      CASE(2)
        e_add_e%E  = A%E  + B%E
        e_add_e%V  = A%V  + B%V
        e_add_e%T2 = A%T2 + B%T2
      CASE(3)
        e_add_e%E  = A%E  + B%E
        e_add_e%V  = A%V  + B%V
        e_add_e%T2 = A%T2 + B%T2
        e_add_e%T3 = A%T3 + B%T3
      CASE(4)
        e_add_e%E  = A%E  + B%E
        e_add_e%V  = A%V  + B%V
        e_add_e%T2 = A%T2 + B%T2
        e_add_e%T3 = A%T3 + B%T3
        e_add_e%T4 = A%T4 + B%T4
      END SELECT 

end function e_add_e

elemental function e_add_dp(A,B)  
   USE PROBLEM_DATA
   type( eb ),intent(in)::  A
   real( dp ),intent(in)::  B
   type( eb )           ::e_add_dp

   select case( deriv_order )
      case(0)
        e_add_dp%E  = A%E  + B
      case(1)
        e_add_dp%E  = A%E  + B
        e_add_dp%V  = A%V 
      CASE(2)
        e_add_dp%E  = A%E  + B
        e_add_dp%V  = A%V  
        e_add_dp%T2 = A%T2 
      CASE(3)
        e_add_dp%E  = A%E  + B
        e_add_dp%V  = A%V  
        e_add_dp%T2 = A%T2 
        e_add_dp%T3 = A%T3 
      CASE(4)
        e_add_dp%E  = A%E  + B
        e_add_dp%V  = A%V  
        e_add_dp%T2 = A%T2 
        e_add_dp%T3 = A%T3 
        e_add_dp%T4 = A%T4 
      END SELECT 

end function e_add_dp

elemental function dp_add_e(A,B)  
   USE PROBLEM_DATA
   real( dp ),intent(in)::  A
   type( eb ),intent(in)::  B
   type( eb )           ::  dp_add_e

   select case( deriv_order )
      case(0)
        dp_add_e%E  = A + B%E  
      case(1)
        dp_add_e%E  = A + B%E  
        dp_add_e%V  =     B%V 
      CASE(2)
        dp_add_e%E  = A + B%E  
        dp_add_e%V  =     B%V  
        dp_add_e%T2 =     B%T2 
      CASE(3)
        dp_add_e%E  = A + B%E  
        dp_add_e%V  =     B%V  
        dp_add_e%T2 =     B%T2 
        dp_add_e%T3 =     B%T3 
      CASE(4)
        dp_add_e%E  = A + B%E  
        dp_add_e%V  =     B%V  
        dp_add_e%T2 =     B%T2 
        dp_add_e%T3 =     B%T3 
        dp_add_e%T4 =     B%T4 
      END SELECT 

end function dp_add_e

!.....DEFINE SUB-OPERATOR RULES FOR (-)

elemental function pv_sub_pv(A,B)
   type( pv ), intent(in )::a,b
   type( pv )             ::pv_sub_pv
   integer::i
   do i = 1, nv
      pv_sub_pv%vpart(I) = a%vpart(I) - b%vpart(I)
   end do
end function pv_sub_pv

elemental function pt2_sub_pt2(A,B)
   type( pt2 ), intent(in )::a,b
   type( pt2 )             ::pt2_sub_pt2
   integer::i
   do i = 1, nv
      pt2_sub_pt2%t2part(:,I) = a%t2part(:,I) - b%t2part(:,I)
   end do
end function pt2_sub_pt2

elemental function pt3_sub_pt3(A,B)
   type( pt3 ), intent(in )::a,b
   type( pt3 )             ::pt3_sub_pt3
   integer::i
   do i = 1, nv
      pt3_sub_pt3%t3part(:,:,I) = a%t3part(:,:,I) - b%t3part(:,:,I)
   end do
end function pt3_sub_pt3

elemental function pt4_sub_pt4(A,B)
   type( pt4 ), intent(in )::a,b
   type( pt4 )             ::pt4_sub_pt4
   integer::i
   do i = 1, nv
      pt4_sub_pt4%t4part(:,:,:,I) = a%t4part(:,:,:,I) - b%t4part(:,:,:,I)
   end do
end function pt4_sub_pt4

!.....DEFINE OCEA OPERATOR RULES FOR (-)

elemental function e_sub_e(A,B)  
   USE PROBLEM_DATA
   type( eb ),intent(in )::A,B
   type( eb )            ::e_sub_e

   select case( deriv_order )
      case(0)
        e_sub_e%E  = A%E - B%E
      case(1)
        e_sub_e%E  = A%E - B%E
        e_sub_e%V  = A%V - B%V
      CASE(2)
        e_sub_e%E  = A%E  - B%E
        e_sub_e%V  = A%V  - B%V
        e_sub_e%T2 = A%T2 - B%T2
      CASE(3)
        e_sub_e%E  = A%E  - B%E
        e_sub_e%V  = A%V  - B%V
        e_sub_e%T2 = A%T2 - B%T2
        e_sub_e%T3 = A%T3 - B%T3
      CASE(4)
        e_sub_e%E  = A%E  - B%E
        e_sub_e%V  = A%V  - B%V
        e_sub_e%T2 = A%T2 - B%T2
        e_sub_e%T3 = A%T3 - B%T3
        e_sub_e%T4 = A%T4 - B%T4
      END SELECT 

end function e_sub_e

!.....DEFINE SUB-OPERATOR RULES FOR (-A)

elemental function minus_pv(A)
   type( pv ), intent(in )::a
   type( pv )             ::minus_pv
   integer::i
   do i = 1, nv
      minus_pv%vpart(I) = -a%vpart(I)
   end do
end function minus_pv

elemental function minus_pt2(A)
   type( pt2 ), intent(in )::a
   type( pt2 )             ::minus_pt2
   integer::i
   do i = 1, nv
      minus_pt2%t2part(:,I) = -a%t2part(:,I)
   end do
end function minus_pt2

elemental function minus_pt3(A)
   type( pt3 ), intent(in )::a
   type( pt3 )             ::minus_pt3
   integer::i
   do i = 1, nv
      minus_pt3%t3part(:,:,I) = -a%t3part(:,:,I)
   end do
end function minus_pt3

elemental function minus_pt4(A)
   type( pt4 ), intent(in )::a
   type( pt4 )             ::minus_pt4
   integer::i
   do i = 1, nv
      minus_pt4%t4part(:,:,:,I) = -a%t4part(:,:,:,I)
   end do
end function minus_pt4

!.....DEFINE OCEA OPERATOR RULES FOR (-)

elemental function minus_e(A)
   USE PROBLEM_DATA
   type( eb ),intent(in )::A
   type( eb )            ::minus_e

   select case( deriv_order )
      case(0)
        minus_e%E  = -A%E
      case(1)
        minus_e%E  = -A%E
        minus_e%V  = -A%V
      CASE(2)
        minus_e%E  = -A%E
        minus_e%V  = -A%V
        minus_e%T2 = -A%T2
      CASE(3)
        minus_e%E  = -A%E
        minus_e%V  = -A%V
        minus_e%T2 = -A%T2
        minus_e%T3 = -A%T3
      CASE(4)
        minus_e%E  = -A%E
        minus_e%V  = -A%V
        minus_e%T2 = -A%T2
        minus_e%T3 = -A%T3
        minus_e%T4 = -A%T4
      END SELECT 

end function minus_e

elemental function e_minus_dp(A,B)
   USE PROBLEM_DATA
   type( eb ),intent(in ):: A
   REAL( dp ),intent(in ):: B
   type( eb )            :: e_minus_dp

   select case( deriv_order )
      case(0)
        e_minus_dp%E  = A%E - B
      case(1)
        e_minus_dp%E  = A%E - B
        e_minus_dp%V  = A%V
      CASE(2)
        e_minus_dp%E  = A%E - B
        e_minus_dp%V  = A%V
        e_minus_dp%T2 = A%T2
      CASE(3)
        e_minus_dp%E  = A%E - B
        e_minus_dp%V  = A%V
        e_minus_dp%T2 = A%T2
        e_minus_dp%T3 = A%T3
      CASE(4)
        e_minus_dp%E  = A%E - B
        e_minus_dp%V  = A%V
        e_minus_dp%T2 = A%T2
        e_minus_dp%T3 = A%T3
        e_minus_dp%T4 = A%T4
      END SELECT 

end function e_minus_dp

elemental function dp_minus_e(A,B)
   USE PROBLEM_DATA
   REAL( dp ),intent(in ):: A
   type( eb ),intent(in ):: B
   type( eb )            :: dp_minus_e

   select case( deriv_order )
      case(0)
        dp_minus_e%E  = A - B%E 
      case(1)
        dp_minus_e%E  = A - B%E 
        dp_minus_e%V  =   - B%V
      CASE(2)
        dp_minus_e%E  = A - B%E 
        dp_minus_e%V  =   - B%V
        dp_minus_e%T2 =   - B%T2
      CASE(3)
        dp_minus_e%E  = A - B%E 
        dp_minus_e%V  =   - B%V
        dp_minus_e%T2 =   - B%T2
        dp_minus_e%T3 =   - B%T3
      CASE(4)
        dp_minus_e%E  = A - B%E 
        dp_minus_e%V  =   - B%V
        dp_minus_e%T2 =   - B%T2
        dp_minus_e%T3 =   - B%T3
        dp_minus_e%T4 =   - B%T4
      END SELECT 

end function dp_minus_e

!.....DEFINE SUB-OPERATOR RULES FOR (S*PT)

elemental function dp_mult_pv(A,B)
   real( dp ), intent(in )::A
   type( pv ), intent(in )::B
   type( pv )             ::dp_mult_pv
   integer::i
   do i = 1, nv
      dp_mult_pv%vpart(I) = a*b%vpart(I)
   end do
end function dp_mult_pv

elemental function dp_mult_pt2(A,B)
   real( dp  ), intent(in )::A
   type( pt2 ), intent(in )::B
   type( pt2 )             ::dp_mult_pt2
   integer::i
   do i = 1, nv
      dp_mult_pt2%t2part(:,I) = a*b%t2part(:,I)
   end do
end function dp_mult_pt2

elemental function dp_mult_pt3(A,B)
   real( dp  ), intent(in )::A
   type( pt3 ), intent(in )::B
   type( pt3 )             ::dp_mult_pt3
   integer::i
   do i = 1, nv
      dp_mult_pt3%t3part(:,:,I) = a*b%t3part(:,:,I)
   end do
end function dp_mult_pt3

elemental function dp_mult_pt4(A,B)
   real( dp  ), intent(in )::A
   type( pt4 ), intent(in )::B
   type( pt4 )             ::dp_mult_pt4
   integer::i
   do i = 1, nv
      dp_mult_pt4%t4part(:,:,:,I) = a*b%t4part(:,:,:,I)
   end do
end function dp_mult_pt4

!.....DEFINE OCEA OPERATOR FOR (S*E)

function DP_mult_E(A,B)
   USE PROBLEM_DATA
   real( dp ),intent(in )::A
   type( eb ),intent(in )::B
   type( eb )            ::dp_mult_e

   select case( deriv_order )
      case(0)
        dp_mult_e%E  = A*B%E
      case(1)
        dp_mult_e%E  = A*B%E
        dp_mult_e%V  = A*B%V
      CASE(2)
        dp_mult_e%E  = A*B%E
        dp_mult_e%V  = A*B%V
        dp_mult_e%T2 = A*B%T2
      CASE(3)
        dp_mult_e%E  = A*B%E
        dp_mult_e%V  = A*B%V
        dp_mult_e%T2 = A*B%T2
        dp_mult_e%T3 = A*B%T3
      CASE(4)
        dp_mult_e%E  = A*B%E
        dp_mult_e%V  = A*B%V
        dp_mult_e%T2 = A*B%T2
        dp_mult_e%T3 = A*B%T3
        dp_mult_e%T4 = A*B%T4
      END SELECT 

end function dp_mult_e

!.....DEFINE SUB-OPERATOR RULES FOR (E*S)

elemental function pv_mult_dp(A,B)
   type( pv ), intent(in )::A
   real( dp ), intent(in )::B
   type( pv )             ::pv_mult_dp
   integer::i
   do i = 1, nv
      pv_mult_dp%vpart(I) = b*a%vpart(I)
   end do
end function pv_mult_dp

elemental function pt2_mult_dp(A,B)
   type( pt2 ), intent(in )::A
   real( dp  ), intent(in )::B
   type( pt2 )             ::pt2_mult_dp
   integer::i
   do i = 1, nv
      pt2_mult_dp%t2part(:,I) = b*a%t2part(:,I)
   end do
end function pt2_mult_dp

elemental function pt3_mult_dp(A,B)
   type( pt3 ), intent(in )::A
   real( dp  ), intent(in )::B
   type( pt3 )             ::pt3_mult_dp
   integer::i
   do i = 1, nv
      pt3_mult_dp%t3part(:,:,I) = b*a%t3part(:,:,I)
   end do
end function pt3_mult_dp

elemental function pt4_mult_dp(A,B)
   type( pt4 ), intent(in )::A
   real( dp  ), intent(in )::B
   type( pt4 )             ::pt4_mult_dp
   integer::i
   do i = 1, nv
      pt4_mult_dp%t4part(:,:,:,I) = b*a%t4part(:,:,:,I)
   end do
end function pt4_mult_dp

!.....DEFINE OCEA OPERATOR FOR (E*S)

function e_mult_dp(A,B)
   USE PROBLEM_DATA
   type( eb ),intent(in )::A
   real( dp ),intent(in )::B
   type( eb )            ::e_mult_dp

   select case( deriv_order )
      case(0)
        e_mult_dp%E  = B*A%E
      case(1)
        e_mult_dp%E  = B*A%E
        e_mult_dp%V  = B*A%V
      CASE(2)
        e_mult_dp%E  = B*A%E
        e_mult_dp%V  = B*A%V
        e_mult_dp%T2 = B*A%T2
      CASE(3)
        e_mult_dp%E  = B*A%E
        e_mult_dp%V  = B*A%V
        e_mult_dp%T2 = B*A%T2
        e_mult_dp%T3 = B*A%T3
      CASE(4)
        e_mult_dp%E  = B*A%E
        e_mult_dp%V  = B*A%V
        e_mult_dp%T2 = B*A%T2
        e_mult_dp%T3 = B*A%T3
        e_mult_dp%T4 = B*A%T4
      END SELECT 

end function e_mult_dp

FUNCTION EB_VEC_MULT_DP_SCA(A,B)
     ! EB VECTOR ARRAY FUNCTION
     TYPE(EB), DIMENSION(:), INTENT(IN):: A
     REAL(DP),               INTENT(IN):: B
     TYPE(EB), DIMENSION(SIZE(A))      :: EB_VEC_MULT_DP_SCA
     INTEGER:: I

     ! EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, SIZE(A)
        EB_VEC_MULT_DP_SCA(I) = A(I)*B
     END DO !I

END FUNCTION EB_VEC_MULT_DP_SCA

FUNCTION DP_SCA_MULT_EB_VEC(A,B)
     ! EB VECTOR ARRAY FUNCTION
     REAL(DP),               INTENT(IN):: A
     TYPE(EB), DIMENSION(:), INTENT(IN):: B
     TYPE(EB), DIMENSION(SIZE(B))      :: DP_SCA_MULT_EB_VEC
     INTEGER:: I

     ! EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, SIZE(B)
        DP_SCA_MULT_EB_VEC(I) = A*B(I)
     END DO !I

END FUNCTION DP_SCA_MULT_EB_VEC

FUNCTION DP_MAT_MULT_EB_VEC(A,B)    RESULT(MxV)
     ! EB VECTOR ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:),INTENT(IN)::  A
     TYPE(EB), DIMENSION(:),  INTENT(IN)::  B
     TYPE(EB), DIMENSION(SIZE(B))       :: MxV 
     INTEGER:: I,J

     !BEGIN DIMENSIONAL ERROR CHECKING

     IF( SIZE(A,2) /= SIZE(B) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_MAT_MULT_EB_VEC INCONSISTENT MATRIX-VECTOR PRODUCT'
        WRITE(6,*)'CURRENT MATRIX DIMS =',SIZE(A,1),SIZE(A,2)
        WRITE(6,*)'CURRENT VECTOR DIM  =',SIZE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX(DP)*VECTOR(EB) MULTIPLICATION

     DO I = 1, SIZE(A,1)
        MxV(I) = A(I,1)*B(1)
        DO J = 2, SIZE(A,2)
           MxV(I) = MxV(I) + A(I,j)*B(j)
        END DO !J
     END DO !I

END FUNCTION DP_MAT_MULT_EB_VEC

FUNCTION EB_MAT_MULT_DP_VEC(A,B)    RESULT(MxV)
     ! EB VECTOR ARRAY FUNCTION
     TYPE(EB), DIMENSION(:,:),INTENT(IN)::  A
     REAL(DP), DIMENSION(:),  INTENT(IN)::  B
     TYPE(EB), DIMENSION(SIZE(B))       :: MxV 
     INTEGER:: I,J

     !BEGIN DIMENSIONAL ERROR CHECKING

     IF( SIZE(A,2) /= SIZE(B) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:EB_MAT_MULT_DP_VEC INCONSISTENT MATRIX-VECTOR PRODUCT'
        WRITE(6,*)'CURRENT MATRIX DIMS =',SIZE(A,1),SIZE(A,2)
        WRITE(6,*)'CURRENT VECTOR DIM  =',SIZE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX(EB)*VECTOR(DP) MULTIPLICATION

     DO I = 1, SIZE(A,1)
        MxV(I) = A(I,1)*B(1)
        DO J = 2, SIZE(A,2)
           MxV(I) = MxV(I) + A(I,j)*B(j)
        END DO !J
     END DO !I

END FUNCTION EB_MAT_MULT_DP_VEC

FUNCTION EB_VEC_MULT_EB_SCA(A,B)
     ! EB VECTOR ARRAY FUNCTION
     TYPE(EB), DIMENSION(:), INTENT(IN):: A
     TYPE(EB),               INTENT(IN):: B
     TYPE(EB), DIMENSION(SIZE(A))      :: EB_VEC_MULT_EB_SCA
     INTEGER::I

     ! EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, SIZE(A)
        EB_VEC_MULT_EB_SCA(I) = A(I)*B
     END DO !I

END FUNCTION EB_VEC_MULT_EB_SCA

FUNCTION EB_SCA_MULT_EB_VEC(A,B)
     ! EB VECTOR ARRAY FUNCTION
     TYPE(EB),               INTENT(IN):: A
     TYPE(EB), DIMENSION(:), INTENT(IN):: B
     TYPE(EB), DIMENSION(SIZE(B))      :: EB_SCA_MULT_EB_VEC
     INTEGER::I

     ! EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, SIZE(B)
        EB_SCA_MULT_EB_VEC(I) = A*B(I)
     END DO !I

END FUNCTION EB_SCA_MULT_EB_VEC

FUNCTION EB_MAT_MULT_DP_SCA(A,B)
     ! EB VECTOR ARRAY FUNCTION
     TYPE(EB), DIMENSION(:,:),     INTENT(IN):: A
     REAL(DP),                     INTENT(IN):: B
     TYPE(EB), DIMENSION(SIZE(A,1),SIZE(A,2)):: EB_MAT_MULT_DP_SCA
     INTEGER:: I,J

     ! EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, SIZE(A,1)
        DO J = 1, SIZE(A,2)
           EB_MAT_MULT_DP_SCA(I,J) = A(I,J)*B
        END DO !J
     END DO !I

END FUNCTION EB_MAT_MULT_DP_SCA

FUNCTION DP_SCA_MULT_EB_MAT(A,B)
     ! EB MATRIX ARRAY FUNCTION
     REAL(DP),                     INTENT(IN):: A
     TYPE(EB), DIMENSION(:,:),     INTENT(IN):: B
     TYPE(EB), DIMENSION(SIZE(B,1),SIZE(B,2)):: DP_SCA_MULT_EB_MAT
     INTEGER:: I,J

     ! EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, SIZE(B,1)
        DO J = 1, SIZE(B,2)
           DP_SCA_MULT_EB_MAT(I,J) = A*B(I,J)
        END DO !J
     END DO !I

END FUNCTION DP_SCA_MULT_EB_MAT

FUNCTION EB_MAT_MULT_EB_SCA(A,B)
     ! EB MATRIX ARRAY FUNCTION
     TYPE(EB), DIMENSION(:,:),     INTENT(IN):: A
     TYPE(EB),                     INTENT(IN):: B
     TYPE(EB), DIMENSION(SIZE(A,1),SIZE(A,2)):: EB_MAT_MULT_EB_SCA
     INTEGER:: I,J

     ! EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, SIZE(A,1)
        DO J = 1, SIZE(A,2)
           EB_MAT_MULT_EB_SCA(I,J) = A(I,J)*B
        END DO !J
     END DO !I

END FUNCTION EB_MAT_MULT_EB_SCA

FUNCTION EB_SCA_MULT_EB_MAT(A,B)
     ! EB MATRIX ARRAY FUNCTION
     TYPE(EB),                     INTENT(IN):: A
     TYPE(EB), DIMENSION(:,:),     INTENT(IN):: B
     TYPE(EB), DIMENSION(SIZE(B,1),SIZE(B,2)):: EB_SCA_MULT_EB_MAT
     INTEGER:: I,J

     ! EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, SIZE(B,1)
        DO J = 1, SIZE(B,2)
           EB_SCA_MULT_EB_MAT(I,J) = A*B(I,J)
        END DO !J
     END DO !I

END FUNCTION EB_SCA_MULT_EB_MAT

FUNCTION EB_MAT_MULT_EB_VEC(A,B)     RESULT(MdotV)
     ! EB VECTOR ARRAY FUNCTION
     TYPE(EB), DIMENSION(:,:),     INTENT(IN):: A
     TYPE(EB), DIMENSION(:),       INTENT(IN):: B
     TYPE(EB), DIMENSION(SIZE(A,1))          :: MdotV
     INTEGER:: I,J

     ! EMBEDDED FUNCTION DEFINITIONS

     DO I = 1, SIZE(A,1)
        MdotV(I) = A(I,1)*B(1)
        DO J = 2, SIZE(A,2)
           MdotV(I) = MdotV(I) + A(I,J)*B(J)
        END DO !J
     END DO !I

END FUNCTION EB_MAT_MULT_EB_VEC


FUNCTION EB_MAT_MULT_EB_MAT(A,B)         RESULT(MxM)
     ! EB MATRIX ARRAY FUNCTION
     TYPE(EB), DIMENSION(:,:), INTENT(IN)    ::  A
     TYPE(EB), DIMENSION(:,:), INTENT(IN)    ::  B
     TYPE(EB), DIMENSION(SIZE(A,1),SIZE(B,2)):: MxM
     INTEGER:: I,J,K,NRA,NCA,NRB,NCB

     NRA=SIZE(A,1); NCA=SIZE(A,2); NRB=SIZE(B,1); NCB=SIZE(B,2)

     IF( NCA /= NRB ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:EB_MAT_MULT_EB_MAT MATRIX-MATRIX PRODUCT NOT DEFINED'
        WRITE(6,*)'CURRENT A MATRIX DIMS =',NRA,NCA
        WRITE(6,*)'CURRENT B MATRIX DIMS =',NRB,NCB
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP 'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX(EB)*MATRIX(EB) MULTIPLICATION

     DO I = 1, NRA
        DO J = 1, NCB
           MxM(I,J) = A(I,1)*B(1,J)
           DO K = 2, NCA
              MxM(I,J) = MxM(I,J) + A(I,K)*B(K,J)
           END DO !K
        END DO !J
     END DO !I

END FUNCTION EB_MAT_MULT_EB_MAT

FUNCTION EB_MAT_MULT_DP_MAT(A,B)         RESULT(MxM)
     ! EB MATRIX ARRAY FUNCTION
     TYPE(EB), DIMENSION(:,:), INTENT(IN)    ::  A
     REAL(DP), DIMENSION(:,:), INTENT(IN)    ::  B
     TYPE(EB), DIMENSION(SIZE(A,1),SIZE(B,2)):: MxM
     INTEGER:: I,J,K,NRA,NCA,NRB,NCB

     NRA=SIZE(A,1); NCA=SIZE(A,2); NRB=SIZE(B,1); NCB=SIZE(B,2)

     IF( NCA /= NRB ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:EB_MAT_MULT_DP_MAT MATRIX-MATRIX PRODUCT NOT DEFINED'
        WRITE(6,*)'CURRENT A MATRIX DIMS =',NRA,NCA
        WRITE(6,*)'CURRENT B MATRIX DIMS =',NRB,NCB
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP 'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX(EB)*MATRIX(DP) MULTIPLICATION

     DO I = 1, NRA
        DO J = 1, NCB
           MxM(I,J) = A(I,1)*B(1,J)
           DO K = 2, NCA
              MxM(I,J) = MxM(I,J) + A(I,K)*B(K,J)
           END DO !K
        END DO !J
     END DO !I

END FUNCTION EB_MAT_MULT_DP_MAT

FUNCTION DP_MAT_MULT_EB_MAT(A,B)         RESULT(MxM)
     ! EB MATRIX ARRAY FUNCTION
     REAL(DP), DIMENSION(:,:), INTENT(IN)    ::  A
     TYPE(EB), DIMENSION(:,:), INTENT(IN)    ::  B
     TYPE(EB), DIMENSION(SIZE(A,1),SIZE(B,2)):: MxM
     INTEGER:: I,J,K,NRA,NCA,NRB,NCB

     NRA=SIZE(A,1); NCA=SIZE(A,2); NRB=SIZE(B,1); NCB=SIZE(B,2)

     IF( NCA /= NRB ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_MAT_MULT_EB_MAT MATRIX-MATRIX PRODUCT NOT DEFINED'
        WRITE(6,*)'CURRENT A MATRIX DIMS =',NRA,NCA
        WRITE(6,*)'CURRENT B MATRIX DIMS =',NRB,NCB
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP 'OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX(DP)*MATRIX(EB) MULTIPLICATION

     DO I = 1, NRA
        DO J = 1, NCB
           MxM(I,J) = A(I,1)*B(1,J)
           DO K = 2, NCA
              MxM(I,J) = MxM(I,J) + A(I,K)*B(K,J)
           END DO !K
        END DO !J
     END DO !I

END FUNCTION DP_MAT_MULT_EB_MAT


!.....DEFINE OCEA OPERATOR FOR DIVISION (/)

elemental function e_DIV_DP(A,B)
   USE PROBLEM_DATA
   type( eb ),intent(in )::A
   real( dp ),intent(in )::B
   real( dp )            ::recip
   type( eb )            ::e_div_dp

   recip = 1.0d0/B

   select case( deriv_order )
      case(0)
        e_div_dp%E  = recip*A%E
      case(1)
        e_div_dp%E  = recip*A%E
        e_div_dp%V  = recip*A%V
      CASE(2)
        e_div_dp%E  = recip*A%E
        e_div_dp%V  = recip*A%V
        e_div_dp%T2 = recip*A%T2
      CASE(3)
        e_div_dp%E  = recip*A%E
        e_div_dp%V  = recip*A%V
        e_div_dp%T2 = recip*A%T2
        e_div_dp%T3 = recip*A%T3
      CASE(4)
        e_div_dp%E  = recip*A%E
        e_div_dp%V  = recip*A%V
        e_div_dp%T2 = recip*A%T2
        e_div_dp%T3 = recip*A%T3
        e_div_dp%T4 = recip*A%T4
      END SELECT 

end function e_div_dp

function dp_DIV_e(A,B)
   USE PROBLEM_DATA
   real( dp ), intent( in ):: A
   type( eb ), intent( in ):: B
   type( eb )              :: dp_DIV_e
   real( dp )                              :: b0, R0
   real( dp ), dimension( nv )             :: b1, R1
   real( dp ), dimension( nv, nv )         :: b2, R2
   real( dp ), dimension( nv, nv, nv )     :: b3, R3
   real( dp ), dimension( nv, nv, nv, nv ) :: b4, R4
   integer :: j,k,r

   R0  = A/B%E

  select case( deriv_order )
      case(0)

        dp_DIV_e%E  = R0

      case(1)
        b0    = b%e
        b1(:) = b%v%vpart(:)

        dp_DIV_e%E  = R0
        dp_DIV_e%V%vpart(:)  = - B1(:)*(R0/B0)

      CASE(2)
        b0      = b%e
        b1(:)   = b%v%vpart(:)
        b2(:,:) = b%t2%t2part(:,:)
        R1(:) = - B1(:)*(R0/B0)

        dp_DIV_e%E  = R0
        dp_DIV_e%V%vpart(:)  = R1(:)
        do j = 1, nv
           dp_DIV_e%T2%t2part(:,j) = &
                    - (B2(:,j)*R0 + B1(:)*R1(j) + R1(:)*B1(j))/B0
        end do !j

      CASE(3)
        b0        = b%e
        b1(:)     = b%v%vpart(:)
        b2(:,:)   = b%t2%t2part(:,:)
        b3(:,:,:) = b%t3%t3part(:,:,:)

!.......NOTE: Ratio do loops can not be mixed for array equations to WORK!!!!!

        R1(:) =  - B1(:)*(R0/B0)

        DO J = 1, NV
           R2(:,J) =  ( - B2(:,j)*R0 - B1(:)*R1(j) - R1(:)*B1(j) )/B0
        END DO !J

        dp_DIV_e%E                      = R0
        dp_DIV_e%V%vpart(:)             = R1(:)
        do j = 1, nv
           dp_DIV_e%T2%t2part(:,j)      = R2(:,J)
           do k = 1, nv
              dp_DIV_e%T3%T3part(:,j,k) = ( - b2(:,j)*r1(k) - b1(:)*r2(j,k)  &
                                            - b2(:,k)*r1(j) - b1(j)*r2(:,k)  &
                                            - b1(k)*r2(:,j) - b2(j,k)*r1(:)  &
                                            - b3(:,j,k)*r0                  )/b0
           end do !k
        end do !j

      CASE(4)
        b0          = b%e
        b1(:)       = b%v%vpart(:)
        b2(:,:)     = b%t2%t2part(:,:)
        b3(:,:,:)   = b%t3%t3part(:,:,:)
        b4(:,:,:,:) = b%t4%t4part(:,:,:,:)

!.......NOTE: Ratio do loops can not be mixed for array equations to WORK!!!!!

        R1(:) =  - B1(:)*(R0/B0)

        DO J = 1, NV
           R2(:,J) =  ( - B2(:,j)*R0 - B1(:)*R1(j) - R1(:)*B1(j) )/B0
        END DO !J

        do j = 1, nv
           do k = 1, nv
              r3(:,j,k) = (           - b2(:,j)*r1(k) - b1(:)*r2(j,k)  &
                                      - b2(:,k)*r1(j) - b1(j)*r2(:,k)  &
                                      - b1(k)*r2(:,j) - b2(j,k)*r1(:)  &
                                      - b3(:,j,k)*r0                  )/b0
           end do !K
        end do !J

        dp_DIV_e%E                      = R0
        dp_DIV_e%V%vpart(:)             = R1(:)
        do j = 1, nv
           dp_DIV_e%T2%t2part(:,j)      = R2(:,J)
           do k = 1, nv
              dp_DIV_e%T3%T3part(:,j,k) = R3(:,J,K)
              do r = 1, nv
                 dp_DIV_e%T4%T4part(:,j,k,r) = &
                         -(B3(:,j,k)*R1(r)  + R3(:,j,k)*B1(r) + B4(:,j,k,r)*R0  &
                         + B2(:,j)*R2(k,r)  + R2(:,j)*B2(k,r) + B3(:,j,r)*R1(k) &
                         + R3(:,j,r)*B1(k)  + B2(:,k)*R2(j,r) + R2(:,k)*B2(j,r) &
                         + B1(:)*R3(j,k,r)  + R1(:)*B3(j,k,r) + B2(:,r)*R2(j,k) &
                         + R2(:,r)*B2(j,k)  + B3(:,k,r)*R1(j) + R3(:,k,r)*B1(j) )/B0
              end do !r
           end do !k
        end do !j
      END SELECT 

end function dp_DIV_e

function e_DIV_e( A,B)
   USE PROBLEM_DATA
   type( eb ), intent( in ):: A,B
   type( eb )              :: e_DIV_e
   real( dp )                              :: a0, b0, R0
   real( dp ), dimension( nv )             :: a1, b1, R1
   real( dp ), dimension( nv, nv )         :: a2, b2, R2
   real( dp ), dimension( nv, nv, nv )     :: a3, b3, R3
   real( dp ), dimension( nv, nv, nv, nv ) :: a4, b4
   integer :: j,k,r,i

   R0  = A%E/B%E

  select case( deriv_order )
      case(0)

        e_DIV_e%E  = R0

      case(1)
        a0    = a%e;           b0    = b%e
        a1(:) = a%v%vpart(:);  b1(:) = b%v%vpart(:)

        e_DIV_e%E  = R0
        e_DIV_e%V%vpart(:)  = ( A1(:) - B1(:)*R0 )/B0

      CASE(2)
        a0      = a%e;              b0      = b%e
        a1(:)   = a%v%vpart(:);     b1(:)   = b%v%vpart(:)
        a2(:,:) = a%t2%t2part(:,:); b2(:,:) = b%t2%t2part(:,:)
        R1(:) = ( A1(:) - B1(:)*R0 )/B0

        e_DIV_e%E  = R0
        e_DIV_e%V%vpart(:)  = R1(:)
        do j = 1, nv
           e_DIV_e%T2%t2part(:,j) = &
                    - (B2(:,j)*R0 + B1(:)*R1(j) + R1(:)*B1(j) - A2(:,j))/B0
        end do !j

      CASE(3)
        a0        = a%e;                b0        = b%e
        a1(:)     = a%v%vpart(:);       b1(:)     = b%v%vpart(:)
        a2(:,:)   = a%t2%t2part(:,:);   b2(:,:)   = b%t2%t2part(:,:)
        a3(:,:,:) = a%t3%t3part(:,:,:); b3(:,:,:) = b%t3%t3part(:,:,:)

!.......NOTE: Ratio do loops can not be mixed for array equations to WORK!!!!!

        R1(:) =  ( A1(:) - B1(:)*R0 )/B0

        DO J = 1, NV
           R2(:,J) =  ( A2(:,J) - B2(:,j)*R0 - B1(:)*R1(j) - R1(:)*B1(j) )/B0
        END DO !J

        e_DIV_e%E                      = R0
        e_DIV_e%V%vpart(:)             = R1(:)
        do j = 1, nv
           e_DIV_e%T2%t2part(:,j)      = R2(:,J)
           do k = 1, nv
              e_DIV_e%T3%T3part(:,j,k) = ( a3(:,j,k) - b2(:,j)*r1(k) - b1(:)*r2(j,k)  &
                                         - b2(:,k)*r1(j) - b1(j)*r2(:,k)  &
                                         - b1(k)*r2(:,j) - b2(j,k)*r1(:)  &
                                         - b3(:,j,k)*r0                  )/b0
           end do !k
        end do !j

      CASE(4)
        a0          = a%e;                  b0          = b%e
        a1(:)       = a%v%vpart(:);         b1(:)       = b%v%vpart(:)
        a2(:,:)     = a%t2%t2part(:,:);     b2(:,:)     = b%t2%t2part(:,:)
        a3(:,:,:)   = a%t3%t3part(:,:,:);   b3(:,:,:)   = b%t3%t3part(:,:,:)
        a4(:,:,:,:) = a%t4%t4part(:,:,:,:); b4(:,:,:,:) = b%t4%t4part(:,:,:,:)

!.......NOTE: Ratio do loops can not be mixed for array equations to WORK!!!!!

        R1(:) =  ( A1(:) - B1(:)*R0 )/B0

        DO J = 1, NV
           R2(:,J) =  ( A2(:,J) - B2(:,j)*R0 - B1(:)*R1(j) - R1(:)*B1(j) )/B0
        END DO !J

        do j = 1, nv
           do k = 1, nv
              r3(:,j,k) = ( a3(:,j,k) - b2(:,j)*r1(k) - b1(:)*r2(j,k)  &
                                      - b2(:,k)*r1(j) - b1(j)*r2(:,k)  &
                                      - b1(k)*r2(:,j) - b2(j,k)*r1(:)  &
                                      - b3(:,j,k)*r0                  )/b0
           end do !K
        end do !J

        e_DIV_e%E                      = R0
        e_DIV_e%V%vpart(:)             = R1(:)
        do j = 1, nv
           e_DIV_e%T2%t2part(:,j)      = R2(:,J)
           do k = 1, nv
              e_DIV_e%T3%T3part(:,j,k) = R3(:,J,K)
              do r = 1, nv
                 e_DIV_e%T4%T4part(:,j,k,r) = &
                         -(B3(:,j,k)*R1(r)  + R3(:,j,k)*B1(r) + B4(:,j,k,r)*R0  &
                         + B2(:,j)*R2(k,r)  + R2(:,j)*B2(k,r) + B3(:,j,r)*R1(k) &
                         + R3(:,j,r)*B1(k)  + B2(:,k)*R2(j,r) + R2(:,k)*B2(j,r) &
                         + B1(:)*R3(j,k,r)  + R1(:)*B3(j,k,r) + B2(:,r)*R2(j,k) &
                         + R2(:,r)*B2(j,k)  + B3(:,k,r)*R1(j) + R3(:,k,r)*B1(j) &
                         - A4(:,j,k,r))/B0
              end do !r
           end do !k
        end do !j

      END SELECT 

end function e_DIV_e

function e_MULT_e( A,B)
   USE PROBLEM_DATA
   type( eb ), intent( in ):: A,B
   type( eb )              :: e_MULT_e
   real( dp )                              :: a0, b0
   real( dp ), dimension( nv )             :: a1, b1
   real( dp ), dimension( nv, nv )         :: a2, b2
   real( dp ), dimension( nv, nv, nv )     :: a3, b3
   real( dp ), dimension( nv, nv, nv, nv ) :: a4, b4
   integer :: j,k,r

  select case( deriv_order )
      case(0)
        a0 = a%e;  b0 = b%e

        e_MULT_e%E  = a0*b0

      case(1)
        a0      = a%e;              b0      = b%e
        a1(:)   = a%v%vpart(:);     b1(:)   = b%v%vpart(:)

        e_MULT_e%E           = a0*b0
        e_MULT_e%V%vpart(:)  = a0*b1(:) + b0*a1(:)

      CASE(2)
        a0      = a%e;              b0      = b%e
        a1(:)   = a%v%vpart(:);     b1(:)   = b%v%vpart(:)
        a2(:,:) = a%t2%t2part(:,:); b2(:,:) = b%t2%t2part(:,:)

        e_MULT_e%E           = a0*b0
        e_MULT_e%V%vpart(:)  = a0*b1(:) + b0*a1(:)
        do j = 1, nv
           e_MULT_e%T2%t2part(:,j) = &
                    a1(:)*b1(j) + b1(:)*a1(j) + a0*b2(:,j) + b0*a2(:,j)
        end do !j

      CASE(3)
        a0        = a%e;                b0        = b%e
        a1(:)     = a%v%vpart(:);       b1(:)     = b%v%vpart(:)
        a2(:,:)   = a%t2%t2part(:,:);   b2(:,:)   = b%t2%t2part(:,:)
        a3(:,:,:) = a%t3%t3part(:,:,:); b3(:,:,:) = b%t3%t3part(:,:,:)

        e_MULT_e%E           = a0*b0
        e_MULT_e%V%vpart(:)  = a0*b1(:) + b0*a1(:)
        do j = 1, nv
           e_MULT_e%T2%t2part(:,j) = &
                    a1(:)*b1(j) + b1(:)*a1(j) + a0*b2(:,j) + b0*a2(:,j)
           do k = 1, nv
              e_MULT_e%T3%T3part(:,j,k) = &
                       a2(:,j)*b1(k) + b2(:,j)*a1(k) + a1(:)*b2(j,k) + &
                       b1(:)*a2(j,k) + a2(:,k)*b1(j) + b2(:,k)*a1(j) + &
                       a0*b3(:,j,k)  + b0*a3(:,j,k)
           end do !k
        end do !j

      CASE(4)
        a0          = a%e;                  b0          = b%e
        a1(:)       = a%v%vpart(:);         b1(:)       = b%v%vpart(:)
        a2(:,:)     = a%t2%t2part(:,:);     b2(:,:)     = b%t2%t2part(:,:)
        a3(:,:,:)   = a%t3%t3part(:,:,:);   b3(:,:,:)   = b%t3%t3part(:,:,:)
        a4(:,:,:,:) = a%t4%t4part(:,:,:,:); b4(:,:,:,:) = b%t4%t4part(:,:,:,:)

        e_MULT_e%E           = a0*b0
        e_MULT_e%V%vpart(:)  = a0*b1(:) + b0*a1(:)
        do j = 1, nv
           e_MULT_e%T2%t2part(:,j) = &
                    a1(:)*b1(j) + b1(:)*a1(j) + a0*b2(:,j) + b0*a2(:,j)
           do k = 1, nv
              e_MULT_e%T3%T3part(:,j,k) = &
                       a2(:,j)*b1(k) + b2(:,j)*a1(k) + a1(:)*b2(j,k) + &
                       b1(:)*a2(j,k) + a2(:,k)*b1(j) + b2(:,k)*a1(j) + &
                       a0*b3(:,j,k)  + b0*a3(:,j,k)
              do r = 1, nv
                 e_MULT_e%T4%T4part(:,j,k,r) = &
                          a3(:,j,k)*b1(r) + b3(:,j,k)*a1(r) + a2(:,j)*b2(k,r) + &
                          b2(:,j)*a2(k,r) + a3(:,j,r)*b1(k) + b3(:,j,r)*a1(k) + &
                          a2(:,k)*b2(j,r) + b2(:,k)*a2(j,r) + a1(:)*b3(j,k,r) + &
                          b1(:)*a3(j,k,r) + a2(:,r)*b2(j,k) + b2(:,r)*a2(j,k) + &
                          a3(:,k,r)*b1(j) + b3(:,k,r)*a1(j) + a0*b4(:,j,k,r)  + &
                          b0*a4(:,j,k,r)
              end do !r
           end do !k
        end do !j

      END SELECT 

end function e_MULT_e


END MODULE EB_VARIABLE_HANDLING
