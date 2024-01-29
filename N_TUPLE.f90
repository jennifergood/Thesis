MODULE N_TUPLE
USE PROBLEM_DATA; USE EB_VARIABLE_HANDLING; IMPLICIT NONE
    
    TYPE TENSOR
        !State AND PARAMETER Transition TENSORS
        REAL(DP)                                 :: T0
             !    SCALAR              
        TYPE(EB),ALLOCATABLE,DIMENSION(:)        :: T1
             !    VECTOR
        REAL(DP),ALLOCATABLE,DIMENSION(:,:)      :: S1
        REAL(DP),ALLOCATABLE,DIMENSION(:,:,:)    :: S2
        REAL(DP),ALLOCATABLE,DIMENSION(:,:,:,:)  :: S3
        REAL(DP),ALLOCATABLE,DIMENSION(:,:,:,:,:)  :: S4
             !    P2 - P3 ARE TENSORS FOR STATE     TRANSITION TENSOR STORAGE

    END TYPE
!.....FULL TENSOR STORAGE MODE FOR DATA    

!.....CREATE N_TUPLE STORAGE FOR EACH VARIABLE, REQUIRED STORAGE SPACE ONLY
!     F:= { T0, T1, S1, S2, S3 }       
!
!     ====================================================
!.....INTRINSIC FUNCTION OVERLOADING INTERFACE ASSIGNMENTS
!     ====================================================
    INTERFACE ASSIGNMENT(=)
        MODULE PROCEDURE TUPLE_EQ_DP, TUPLE_EQ_TUPLE,&
                         DP_T2_EQ_TUPLE_SCA, DP_T3_EQ_TUPLE_SCA , DP_T5_EQ_TUPLE_SCA, &
                         TUPLE_SCA_EQ_EB_VEC                     
    END INTERFACE
    
    INTERFACE OPERATOR(+)
        MODULE PROCEDURE TUPLE_ADD_TUPLE
    END INTERFACE
    
    INTERFACE OPERATOR(-)
        MODULE PROCEDURE TUPLE_SUB_TUPLE
    END INTERFACE


    INTERFACE OPERATOR(*)
        MODULE PROCEDURE DP_MULT_NTUPLE
    END INTERFACE
    
    
    !.....DEFINE OVERLOADED EMBEDDED OBJECT OPERATOR FUNCTION INTERFACE

   INTERFACE OPERATOR(.D.)   !CONTRACTION OPERATIONS
     MODULE PROCEDURE                                  &
        EB_VEC_DOT_EB_VEC, EB_MAT_DOT_EB_VEC,          &
        DP_VEC_DOT_DP_VEC, DP_MAT_DOT_DP_VEC,          &
        DP_T3_DOT_DP_VEC,  DP_T4_DOT_DP_VEC,           &
        DP_T5_DOT_DP_VEC,  DP_T2_DOT_DP_T3,            &
        DP_T2_DOT_DP_T2,   DP_T3_DOT_DP_T2,            &
        DP_T4_DOT_DP_T2,   DP_T5_DOT_DP_T2,            &
        DP_T3_DOT_DP_T3,   DP_T3_DOT_DP_T4,            &
        DP_T2_DOT_DP_T4,   DP_T4_DOT_DP_T3,            &
        DP_T2_DOT_DP_T5
   END INTERFACE
   

CONTAINS

   SUBROUTINE ALLOCATE_NTUPLE( A )
    !.....THE SCALAR N-TUPLE DATA STRUCTURE IS ASSUMED TO BE
    !.....GIVEN BY:
    !                 G:= { T,  X,  S1, S2, S3, S4 } 
    !.....X     DENOTES STATE VARIABLES FOR DXDT = F(X,T) (OCEA Nx1)
    !....P2     FIRST  ORDER STATE TRANSITION MATRIX FOR DS1/DT = DEL(F)*S1; S1(T0) = I (DP NxN)
    !....P3     SECOND ORDER STATE TRANSITION MATRIX FOR DS2/DT = DEL(DEL(F))*S1*S1 + DEL(F)*S2; P2(T0) = 0 (DP NxNxN)

    TYPE(TENSOR),INTENT(INOUT):: A
          
        ALLOCATE( A%T1(NV) ); A%T1 = ZERO
        SELECT CASE( P_ORDER )
            CASE( 0 )
                WRITE(*,*); WRITE(*,*)'SETUP: NO STATE TRANSITION TENSOR CALCULATION WILL BE PREFORMED!!!'
            CASE( 1 )
                ALLOCATE( A%S1( NV, NV ) ); A%S1 = ZERO
            CASE( 2 )
                ALLOCATE( A%S1( NV, NV     ) ); A%S1 = ZERO
                ALLOCATE( A%S2( NV, NV, NV ) ); A%S2 = ZERO
            CASE( 3 )
                ALLOCATE( A%S1( NV, NV         ) ); A%S1 = ZERO 
                ALLOCATE( A%S2( NV, NV, NV     ) ); A%S2 = ZERO 
                ALLOCATE( A%S3( NV, NV, NV, NV ) ); A%S3 = ZERO 
            CASE( 4 )
                ALLOCATE( A%S1( NV, NV         ) ); A%S1 = ZERO 
                ALLOCATE( A%S2( NV, NV, NV     ) ); A%S2 = ZERO 
                ALLOCATE( A%S3( NV, NV, NV, NV ) ); A%S3 = ZERO 
                ALLOCATE( A%S4( NV, NV, NV, NV, NV ) ); A%S4 = ZERO 
            CASE DEFAULT
                WRITE(*,*            )'ERROR MESSAGE PRINTED TO ERROR_HISTORY'
                WRITE(ERROR_HISTORY,*)
                WRITE(ERROR_HISTORY,*)'FAILED IN SETUP ALLOCATION FOR P_STM:'
                WRITE(ERROR_HISTORY,*)'    P_ORDER NOT DEFINED FOR STATE TRANSITION TENSOR'
                WRITE(ERROR_HISTORY,*)'P_ORDER =',P_ORDER
                STOP                  'STATE TRANSITION TENSOR ALLOCATION ERROR: SETUP'
        END SELECT

            
    END SUBROUTINE ALLOCATE_NTUPLE

!.....DEFINE SUB-OBJECT ASSIGNMENT RULES FOR (=)

    SUBROUTINE TUPLE_EQ_DP( A, B )
        TYPE(TENSOR),INTENT(INOUT):: A
        REAL(DP    ),INTENT(IN   ):: B
        
        A%T0    = B
        A%T1(:) = B
        
        IF( P_ORDER >=1 ) A%S1(:,:    ) = B
        IF( P_ORDER >=2 ) A%S2(:,:,:  ) = B
        IF( P_ORDER >=3 ) A%S3(:,:,:,:) = B
        IF( P_ORDER >=4 ) A%S4(:,:,:,:,:) = B

    END SUBROUTINE TUPLE_EQ_DP
!
!.....
!
    SUBROUTINE tuple_eq_tuple( A, b ) 
        TYPE(TENSOR),INTENT(INOUT):: A
        type(TENSOR),intent(in   ):: b
        integer                   :: i
                    
        A%T0 = B%T0
        A%T1 = B%T1
        
        IF( P_ORDER >=1 ) A%S1 = B%S1
        IF( P_ORDER >=2 ) A%S2 = B%S2
        IF( P_ORDER >=3 ) A%S3 = B%S3
        IF( P_ORDER >=4 ) A%S4 = B%S4
        
    END SUBROUTINE tuple_eq_tuple
!
!.....
!   
    FUNCTION DP_mult_NTUPLE(A,B) RESULT( DPxTUPLE )
        real(  dp  ),intent(in):: A
        TYPE(TENSOR),INTENT(IN):: B
        TYPE(TENSOR)           :: DPxTUPLE
        
        CALL ALLOCATE_NTUPLE( DPxTUPLE ); DPxTUPLE = ZERO
                  
        DPxTUPLE%T0 = A*B%T0
        DPxTUPLE%T1 = A*B%T1
        
        IF( P_ORDER >=1 ) DPxTUPLE%S1 = A*B%S1
        IF( P_ORDER >=2 ) DPxTUPLE%S2 = A*B%S2
        IF( P_ORDER >=3 ) DPxTUPLE%S3 = A*B%S3
        IF( P_ORDER >=4 ) DPxTUPLE%S4 = A*B%S4

    end FUNCTION DP_MULT_NTUPLE
!
!.....
!
    FUNCTION TUPLE_ADD_TUPLE( a, b ) RESULT( AaddB )
        use problem_data 
        type(TENSOR),intent(IN ):: A,B
        type(TENSOR)            :: AaddB
        LOGICAL(lgt)            :: FPASS = .TRUE.
        
        CALL ALLOCATE_NTUPLE( AaddB ); AaddB = ZERO
               
        AaddB%T0 = A%T0 + B%T0
        AaddB%T1 = A%T1 + B%T1
        
        IF( P_ORDER >=1 ) AaddB%S1 = A%S1 + B%S1
        IF( P_ORDER >=2 ) AaddB%S2 = A%S2 + B%S2
        IF( P_ORDER >=3 ) AaddB%S3 = A%S3 + B%S3
        IF( P_ORDER >=4 ) AaddB%S4 = A%S4 + B%S4
                
    END FUNCTION TUPLE_ADD_TUPLE
!
!.....
!
    FUNCTION TUPLE_SUB_TUPLE( a, b ) RESULT( AsubB )
        use problem_data 
        type(TENSOR),intent(IN ):: A,B
        type(TENSOR)            :: AsubB

        CALL ALLOCATE_NTUPLE( AsubB ); AsubB = ZERO
            
        AsubB%T0 = A%T0 - B%T0
        AsubB%T1 = A%T1 - B%T1
        
        IF( P_ORDER >=1 ) AsubB%S1 = A%S1 - B%S1
        IF( P_ORDER >=2 ) AsubB%S2 = A%S2 - B%S2
        IF( P_ORDER >=3 ) AsubB%S3 = A%S3 - B%S3
        IF( P_ORDER >=4 ) AsubB%S4 = A%S4 - B%S4
        
    END FUNCTION TUPLE_SUB_TUPLE
!
!.....
!  
    SUBROUTINE TUPLE_SCA_EQ_EB_VEC(A,B)
        TYPE(TENSOR),             INTENT(INOUT):: A
        TYPE(EB    ),DIMENSION(:),INTENT(IN   ):: B
        
        A%T1(:) = B
        
    END SUBROUTINE TUPLE_SCA_EQ_EB_VEC
    
    SUBROUTINE DP_T2_EQ_TUPLE_SCA(A,B)  
        REAL(DP    ),DIMENSION(:,:),INTENT(OUT):: A
        TYPE(TENSOR),               INTENT(IN ):: B

        A(:,:) = B%S1(:,:)  
                           
    END SUBROUTINE DP_T2_EQ_TUPLE_SCA
!
!.....
!    
    SUBROUTINE DP_T3_EQ_TUPLE_SCA(A,B)
        REAL(DP    ),DIMENSION(:,:,:),INTENT(OUT):: A
        TYPE(TENSOR),                 INTENT(IN ):: B
        
        A(:,:,:) = B%S2(:,:,:)
    
    END SUBROUTINE DP_T3_EQ_TUPLE_SCA
!
!.....
!    
    SUBROUTINE DP_T4_EQ_TUPLE_SCA(A,B)
        REAL(DP    ),DIMENSION(:,:,:,:),INTENT(OUT):: A
        TYPE(TENSOR),                   INTENT(IN ):: B
        
        A(:,:,:,:) = B%S3(:,:,:,:)
    
    END SUBROUTINE DP_T4_EQ_TUPLE_SCA
    
    SUBROUTINE DP_T5_EQ_TUPLE_SCA(A,B)
        REAL(DP    ),DIMENSION(:,:,:,:,:),INTENT(OUT):: A
        TYPE(TENSOR),                   INTENT(IN ):: B
        
        A(:,:,:,:,:) = B%S4(:,:,:,:,:)
    
    END SUBROUTINE DP_T5_EQ_TUPLE_SCA
!
!.....
!     
    !...UTILITIES FOR OVERLOADING THE ".DOT.= DOT PRODUCT" ASSIGNMENT OPERATOR

FUNCTION  DP_VEC_DOT_DP_VEC(A,B)  RESULT(VdotV)
     REAL(DP)::                          VdotV
     REAL(DP), DIMENSION(:), INTENT(IN):: A,B
     ! EMBEDDED FUNCTION DEFINITIONS
     VdotV = DOT_PRODUCT( A, B )
END FUNCTION DP_VEC_DOT_DP_VEC

FUNCTION DP_MAT_DOT_DP_VEC(A,B)      RESULT(MdotV)
     REAL(DP), DIMENSION(:,:), INTENT(IN):: A
     REAL(DP), DIMENSION(:  ), INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(B))        :: MdotV
     INTEGER:: I,J
     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-VECTOR MULTIPLICATION

     MdotV = MATMUL( A, B )

END FUNCTION DP_MAT_DOT_DP_VEC

FUNCTION DP_T3_DOT_DP_VEC(A,B)           RESULT(T3dotV)
     ! DP MATRIX ARRAY FUNCTION: Tij = Tijp.Vp
     REAL(DP), DIMENSION(:,:,:),   INTENT(IN):: A
     REAL(DP), DIMENSION(:    ),   INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)):: T3dotV
     INTEGER:: I,J,K

     IF( SIZE(A,3) /= SIZE(B) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T3_DOT_DP_VEC INCONSISTENT MATRIX-VECTOR DIMS'
        WRITE(6,*)'CURRENT MATRIX DIMS =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM  =',SIZE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        WRITE(6,*)'EXPECTED RESULT: SIZE(A,3) = SIZE(B)'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-VECTOR MULTIPLICATION
     
     DO I = 1, SIZE(A,1)
        DO J = 1, SIZE(A,2)
           T3dotV(I,J) = DOT_PRODUCT( A(I,J,:), B(:) )
        END DO !J
     END DO !I

END FUNCTION DP_T3_DOT_DP_VEC

FUNCTION DP_T4_DOT_DP_VEC(A,B)     RESULT( T4dotV )
     ! DP MATRIX ARRAY FUNCTION: Tijk = Tijkp.Vp
     REAL(DP), DIMENSION(:,:,:,:),           INTENT(IN):: A
     REAL(DP), DIMENSION(:    ),             INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(A,3)):: T4dotV
 
     INTEGER::I,J,K,L

     IF( SIZE(A,4) /= SIZE(B) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T4_DOT_DP_VEC INCONSISTENT MATRIX-VECTOR DIMS'
        WRITE(6,*)'CURRENT MATRIX DIMS =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM  =',SIZE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        WRITE(6,*)'EXPECTED RESULT: SIZE(A,4) = SIZE(B)'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
     END IF
     
     DO I = 1, SIZE(A,1)
        DO J = 1, SIZE(A,2)
           DO K = 1, SIZE(A,3)
              T4dotV(I,J,K) = DOT_PRODUCT( A(I,J,K, :), B(:) ) 
           END DO !K
        END DO ! J
     END DO ! I

END FUNCTION DP_T4_DOT_DP_VEC

FUNCTION DP_T5_DOT_DP_VEC(A,B)     RESULT( T5dotV )
     ! DP MATRIX ARRAY FUNCTION: Tijkl = Tijklp.Vp
     REAL(DP), DIMENSION(:,:,:,:,:),                   INTENT(IN):: A
     REAL(DP), DIMENSION(:    ),                       INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(A,3),SIZE(A,4)):: T5dotV
 
     INTEGER::I,J,K,L

     IF( SIZE(A,5) /= SIZE(B) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T5_DOT_DP_VEC INCONSISTENT MATRIX-VECTOR DIMS'
        WRITE(6,*)'CURRENT MATRIX DIMS =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM  =',SIZE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        WRITE(6,*)'EXPECTED RESULT: SIZE(A,5) = SIZE(B)'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
     END IF
     
     DO I = 1, SIZE(A,1)
        DO J = 1, SIZE(A,2)
           DO K = 1, SIZE(A,3)
              DO L = 1, SIZE(A,4)
                 T5dotV(I,J,K,L) = DOT_PRODUCT( A(I,J,K,L,:), B(:) )
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
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     VdotV = A(1)*B(1)
     DO I = 2, SIZE(B)
        VdotV = VdotV + A(I)*B(I)
     END DO

END FUNCTION EB_VEC_DOT_EB_VEC

FUNCTION EB_MAT_DOT_EB_VEC(A,B)      RESULT(MdotV)
     ! EB VECTOR ARRAY FUNCTION
     TYPE(EB), DIMENSION(:,:), INTENT(IN):: A
     TYPE(EB), DIMENSION(:  ), INTENT(IN):: B
     TYPE(EB), DIMENSION(SIZE(A))        :: MdotV
     INTEGER:: I,J

     IF( SIZE(A,2) /= SIZE(B) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:EB_MAT_DOT_EB_VEC DIM A(*,2) /= DIM B'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SIZE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
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
     ! EB VECTOR ARRAY FUNCTION: Tij = Tip.Tpj
     REAL(DP), DIMENSION(:,:),      INTENT(IN):: A
     REAL(DP), DIMENSION(:,:),      INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(B,2)) :: T2dotT2
     INTEGER:: n1,n2

     IF( SIZE(A,2) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T2_DOT_DP_T2 DIM A(*,2) /= DIM B(1,*)'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SHAPE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
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
     ! EB VECTOR ARRAY FUNCTION: Tijk = Tip.Tpjk
     REAL(DP), DIMENSION(:,:  ),             INTENT(IN):: A
     REAL(DP), DIMENSION(:,:,:),             INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(B,2),SIZE(B,3)):: T2dotT3
     INTEGER:: n1,n2,n3

     IF( SIZE(A,2) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T2_DOT_DP_T3 DIM A(*,2) /= DIM B(1,*,*)'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SHAPE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
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
     ! EB VECTOR ARRAY FUNCTION: Tijk = Tijp.Tpk
     REAL(DP), DIMENSION(:,:,:),             INTENT(IN):: A
     REAL(DP), DIMENSION(:,:  ),             INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(B,2)):: T3dotT2
     INTEGER:: n1,n2,n3

     IF( SIZE(A,3) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T3_DOT_DP_T2 DIM A(*,*,3) /= DIM B(1,*)'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SHAPE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
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
     ! EB VECTOR ARRAY FUNCTION: Tijkl = Tijkp.Tpl
     REAL(DP), DIMENSION(:,:,:,:),                     INTENT(IN):: A
     REAL(DP), DIMENSION(:,:    ),                     INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(A,3),SIZE(B,2)):: T4dotT2
     INTEGER:: n1,n2,n3,n4

     IF( SIZE(A,4) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T4_DOT_DP_T2 DIM A(*,*,*,4) /= DIM B(1,*)'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SHAPE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
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

FUNCTION DP_T5_DOT_DP_T2(A,B)      RESULT(T5dotT2)
     ! EB VECTOR ARRAY FUNCTION: Tijkl = Tijkp.Tpl
     REAL(DP), DIMENSION(:,:,:,:,:),                   INTENT(IN):: A
     REAL(DP), DIMENSION(:,:    ),                     INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(A,3),SIZE(A,4),SIZE(B,2)):: T5dotT2
     INTEGER:: n1,n2,n3,n4,n5

     IF( SIZE(A,5) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T5_DOT_DP_T2 DIM A(*,*,*,4) /= DIM B(1,*)'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SHAPE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-TENSOR MULTIPLICATION A(I,J,K,:)*B(:,M)
     
     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(A,2)
           DO n3=1,SIZE(A,3)
              DO n4=1,SIZE(A,4) 
                 DO n5=1,SIZE(B,2)
                    T5dotT2(n1,n2,n3,n4,n5) = dot_product( A(n1,n2,n3,n4,:), B(:,n5) )
                 END DO   
              END DO
           END DO
        END DO
     END DO

END FUNCTION DP_T5_DOT_DP_T2

FUNCTION DP_T2_DOT_DP_T4(A,B)      RESULT(T2dotT4)
     ! EB VECTOR ARRAY FUNCTION: Tijkl = Tip.Tpjkl
     REAL(DP), DIMENSION(:,:    ),                     INTENT(IN):: A
     REAL(DP), DIMENSION(:,:,:,:),                     INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(B,2),SIZE(B,3),SIZE(B,4)):: T2dotT4
     INTEGER:: n1,n2,n3,n4

     IF( SIZE(A,2) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T2_DOT_DP_T4 DIM A(*,2) /= DIM B(1,*,*,*)'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SHAPE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
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

FUNCTION DP_T2_DOT_DP_T5(A,B)      RESULT(T2dotT5)
     ! EB VECTOR ARRAY FUNCTION: Tijkl = Tip.Tpjkl
     REAL(DP), DIMENSION(:,:    ),                     INTENT(IN):: A
     REAL(DP), DIMENSION(:,:,:,:,:),                   INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(B,2),SIZE(B,3),SIZE(B,4),SIZE(B,5)):: T2dotT5
     INTEGER:: n1,n2,n3,n4,n5

     IF( SIZE(A,2) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T2_DOT_DP_T5 DIM A(*,2) /= DIM B(1,*,*,*)'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SHAPE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     !...MATRIX-TENSOR MULTIPLICATION A(I,:)*B(:,K,L,M)
     
     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(B,2)
           DO n3=1,SIZE(B,3)
              DO n4=1,SIZE(B,4)
                 DO n5=1,SIZE(B,5)
                    T2dotT5(n1,n2,n3,n4,n5) = dot_product( A(n1,:), B(:,n2,n3,n4,n5) )
                 END DO
              END DO
           END DO
        END DO
     END DO

END FUNCTION DP_T2_DOT_DP_T5

FUNCTION DP_T3_DOT_DP_T3(A,B)      RESULT(T3dotT3)
     ! EB VECTOR ARRAY FUNCTION: Tijkl = Tijp.Tpkl
     REAL(DP), DIMENSION(:,:,:),                       INTENT(IN):: A
     REAL(DP), DIMENSION(:,:,:),                       INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(B,2),SIZE(B,3)):: T3dotT3
     INTEGER:: n1,n2,n3,n4

     IF( SIZE(A,3) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T3_DOT_DP_T3 DIM A(*,*,3) /= DIM B(1,*,*)'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SHAPE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(A,2)
           DO n3=1,SIZE(B,2)
              DO n4=1,SIZE(B,3)
                 T3dotT3(n1,n2,n3,n4) = dot_product( A(n1,n2,:), B(:,n3,n4) )
              END DO
           END DO
        END DO
     END DO

END FUNCTION DP_T3_DOT_DP_T3 

FUNCTION DP_T3_DOT_DP_T4(A,B)      RESULT(T3dotT4)
     ! EB VECTOR ARRAY FUNCTION: Tijkl = Tijp.Tpkl
     REAL(DP), DIMENSION(:,:,:),                       INTENT(IN):: A
     REAL(DP), DIMENSION(:,:,:,:),                     INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(B,2),SIZE(B,3),SIZE(B,4)):: T3dotT4
     INTEGER:: n1,n2,n3,n4,n5

     IF( SIZE(A,3) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T3_DOT_DP_T4 DIM A(*,*,3) /= DIM B(1,*,*)'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SHAPE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(A,2)
           DO n3=1,SIZE(B,2)
              DO n4=1,SIZE(B,3)
                  DO n5=1,SIZE(B,4)
                     T3dotT4(n1,n2,n3,n4,n5) = dot_product( A(n1,n2,:), B(:,n3,n4,n5) )
                  END DO
              END DO
           END DO
        END DO
     END DO

END FUNCTION DP_T3_DOT_DP_T4 

FUNCTION DP_T4_DOT_DP_T3(A,B)      RESULT(T4dotT3)
     ! EB VECTOR ARRAY FUNCTION: Tijkl = Tijp.Tpkl
     REAL(DP), DIMENSION(:,:,:,:),                     INTENT(IN):: A
     REAL(DP), DIMENSION(:,:,:),                       INTENT(IN):: B
     REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2),SIZE(A,3),SIZE(B,2),SIZE(B,3)):: T4dotT3
     INTEGER:: n1,n2,n3,n4,n5

     IF( SIZE(A,4) /= SIZE(B,1) ) THEN
        WRITE(6,*)
        WRITE(6,*)'ERROR:DP_T4_DOT_DP_T3 DIM A(*,*,3) /= DIM B(1,*,*)'
        WRITE(6,*)'CURRENT MATRIX DIM =',SHAPE(A)
        WRITE(6,*)'CURRENT VECTOR DIM =',SHAPE(B)
        WRITE(6,*)'FIX PREVIOUSLY IDENTIFIED DIMENSION ERRORS'
        STOP      'N_TUPLE: OPERATOR OVERLOADING ERRORS DETECTED'
     END IF

     ! EMBEDDED FUNCTION DEFINITIONS

     DO n1=1,SIZE(A,1)
        DO n2=1,SIZE(A,2)
           DO n3=1,SIZE(A,2)
              DO n4=1,SIZE(B,2)
                 DO n5=1,SIZE(B,3)
                    T4dotT3(n1,n2,n3,n4,n5) = dot_product( A(n1,n2,n3,:), B(:,n4,n5) )
                 END DO
              END DO
           END DO
        END DO
     END DO

END FUNCTION DP_T4_DOT_DP_T3 
 
END MODULE N_TUPLE