
    
    
        MODULE INTEGRATION
USE PROBLEM_DATA
USE N_TUPLE
USE EB_VARIABLE_HANDLING
USE OCEA_IO
USE OCEA_MATRIX_LIB
IMPLICIT NONE

! SUPPORTED ROUTINES FOR NUMERICAL INTEGRATION
!=====================================================================
! SETUP                     EB,      GENERATE OCEA IC'S
! INIT_COND                 EB,      DEFINE INITIAL CONDITIONS
! EB_RUNGE_KUTTA_DRIVER     EB,      RK4 DRIVER FOR INTEGRATION
! EB_RK4                    EB,      RK4 DERIVATIVE CALLS
! DERIVATIVE		        EB,      VECTOR DERIVATIVES FOR INTEGRATION
! fix_d0_dC                 DP,      COMPUTE INTEGRATION WEIGHTS
! LAST_STEP		       LOGICAL,      TEST END OF INTEGRATION & RESET DT

!.....THIS PROGRAM INTEGRATES A GENERALIZED NTUPLE SCALAR VARIABLE.  OCEA IS 
!     USED TO GENERATE THE PARTIAL DERIVATIVES REQUIRED FOR COMPUTING STATE
!     AND PARAMETER STATE TRANSITION TENSORS.  
CONTAINS
!
!=================================================================
!
SUBROUTINE INIT_COND( X0 )
    !*****************************************************************
    !
    !..PURPOSE:     THIS PROGRAM COMPUTES THE INITIAL CONDITIONS
    !               FOR THE NUMERICAL INTEGRATION PROCESS
    !=======================================================================
    !COPYRIGHT(C) 2009 JAMES D. TURNER, ALL RIGHTS RESERVED
    !=======================================================================
    IMPLICIT NONE
    
    TYPE(EB    ), DIMENSION(:)  , INTENT(INOUT):: X0
              !OCEA EMBEDDED STATE VECTOR, DOF = NV+NP 
                                                        !1234567890123456789012345678901
    CHARACTER(LEN=31)                          :: PATH = 'MODULE_INTEGRATION,INIT_COND:'              
                          
    INTEGER :: I
  


       !.....APPLICATION PARAMETER STATEMENTS
    real(dp    ),parameter:: X_e     = 6780.0d3  ! KM
    real(dp    ),parameter:: Y_e     = 100.0d3   ! KM
    real(dp    ),parameter:: Z_e     = 200.0d3   ! KM
    real(dp    ),parameter:: dX_e    = -0.1d3     !inertial velocity KM/S
    real(dp    ),parameter:: dY_e    =  8.1d3     !inertial velocity KM/S
    real(dp    ),parameter:: dZ_e    = -0.2d3     !inertial velocity KM/S


    
    CHECK_MAX_DERIVATIVE_ORDER: IF( DERIV_ORDER>=5 ) THEN 
         WRITE(*,*) 
         WRITE(*,*)'TENSOR_RATES: ERROR WRITE_TRANSITION_TENSOR_CODE'
         WRITE(*,*)'ERROR: EXP_ORDER>=5 =>INCREASE EXP_ORDER IN PROBLEM DATA'
         STOP 'TENSOR'
    END IF CHECK_MAX_DERIVATIVE_ORDER
            
!.....INTEGRATION START AND STOP TIMES( T0 & TF): DEFINED IN PROBLEM_DATA
    
!.....NUMBER OF INTEGRATION STEPS (NSTEPS):       DEFINED IN PROBLEM_DATA
    
      WRITE(*,*); WRITE(*,rFMT)PATH//' DT =',DT 
!.....ESTIMATED INTEGRATION STEP SIZE (DT):       DEFINED IN PROBLEM_DATA

!.....DEFINE INITIAL CONDITIONS           
!.....NEED TO READ IN STATE INITIAL CONDITION FILE!!!!!!!!!
      X0 = .CONSTANT.ZERO
      DO I = 1, NV; X0(I)%V%VPART(I) = ONE; END DO
      
!.....INITIAL STATE VECTOR AND DERIVATIVES

           X0(1)%E = X_e
           X0(2)%E = Y_e
           X0(3)%E = Z_e
           X0(4)%E = dX_e 
           X0(5)%E = dY_e
           X0(6)%E = dZ_e
!              LOAD INITIAL STATE VECTOR

           write(TIME_HISTORY,*)"     time       x0(1)       x0(2)  ...."
           write(TIME_HISTORY,'(F12.5,10E20.11)')time,(X0(i)%E,i=1,nv)
               
           
           write(StateTime,*)"     time       x0(1)       x0(2)  ...."
           write(StateTime,'(F12.5,10E20.11)')time,(X0(i)%E,i=1,nv)
           
!NOTE: STATE VECTOR AND OCEA INDEPENDENT VARIABLES ELEMENTS ARE LOADED ELSEWHERE

    WRITE(ECHO_INPUT,*)
    WRITE(ECHO_INPUT,iFMT) PATH
    WRITE(ECHO_INPUT,iFMT)'NSTEPS =',NSTEPS
    WRITE(ECHO_INPUT,rFMT)'T0,TF  =',T0,TF
    WRITE(ECHO_INPUT,rFMT)'DT     =',DT
    WRITE(ECHO_INPUT,iFMT)'NV     =',NV
    
    WRITE(*,*); WRITE(*,*) PATH
    WRITE(*,*)'OUTPUT PRINTED TO FILE: TIME_HISTORY.DAT (FILE=30)'
    WRITE(*,*) PATH//'FINISHED INN INIT_COND'

END SUBROUTINE INIT_COND
!
!=================================================================
!
SUBROUTINE RUNGE_KUTTA_DRIVER( X )
   USE problem_data
   USE OCEA_IO; USE OCEA_MATRIX_LIB
    IMPLICIT NONE
    !*****************************************************************
    !
    !..PURPOSE:     THIS PROGRAM INTEGRATES THE VECTOR SET OF
    !               EQUATION (X0) FROM { T0 => TF }.  THE CLASSICAL
    !               4-TH ORDER SET OF EQUATIONS ARE USED.
    !*****************************************************************
    !.....ARGUMENT LIST VARIABLES
    !
    !__CALLED BY: MAINPROGRAM
    !*****************************************************************
    !.....THE SCALAR N-TUPLE DATA STRUCTURE IS ASSUMED TO BE
    !.....GIVEN BY:
    !                 G:= { T,  X,  P2, P3, P4 } 
    !
    !.....T     DENOTES TIME IN THE SIMULATION
    !
    !.....X     DENOTES STATE VARIABLES FOR 
    !                DXDT   = F(X,T) (OCEA Nx1)
    !
    !....P2     FIRST  ORDER STATE TRANSITION MATRIX FOR 
    !                DP2/DT = DEL(F)*P2; P2(T0) = I (DP NxN)
    !
    !....P3     SECOND ORDER STATE TRANSITION MATRIX FOR 
    !                DP3/DT = DEL(DEL(F))P2*P2+DEL(F)*P3; P3(T0) = 0 (DP NxNxN)
    !
    !....P4     THIRD ORDER STATE TRANSITION MATRIX FOR 
    !                DP4/DT = DEL(DEL(DEL(F)))*P2*P2*P2 + 2*DEL(DEL(F))*P3*P2 + &
    !                DEL(DEL(F)*P2*P3 + DEL(F)*P4       ; P4(T0) = 0 (DP NxNxNxN)
    !
    !=======================================================================
    !COPYRIGHT(C) 2009 JAMES D. TURNER, ALL RIGHTS RESERVED
    !AMENDED: 2011
    !=======================================================================

   
    
    logical                  :: done      
                               !INTEGRATION LOGIC VARIABLE FOR FINAL INTEGRATION STEP
    !LOGICAL,  PARAMETER:: P_STM       = .TRUE.
                                            
    INTEGER                  :: STEP      
                               !NUMBER OF INTEGRATION STEPS TAKEN
                                            
    REAL(DP)                 :: DELT, SUM      
                               !INTEGRATION STEP SIZE

    REAL(DP), ALLOCATABLE,&
              DIMENSION(:,:) :: DXDX0
     
    TYPE(EB    ), DIMENSION(:) , INTENT(INOUT):: X
      
    TYPE(EB    ),DIMENSION(NV):: xc                          

        TYPE(EB    ),DIMENSION(NV):: xp                          

        TYPE(EB    ),DIMENSION(NV):: dxdt                          

                               !FIRST ORDER: STATE TRANSITION MATRICES  
                               
    REAL(DP), DIMENSION(NV,NV,NV,NV,NV) :: Sumplectic_4th_order_STM
    
    REAL(DP), DIMENSION(NV,NV) ::STM_1ST
    REAL(DP), DIMENSION(NV,NV,NV) ::STM_2ND
    REAL(DP), DIMENSION(NV,NV,NV,NV) ::STM_3RD
    REAL(DP), DIMENSION(NV,NV,NV,NV,NV) ::STM_4TH
    
    REAL(DP), DIMENSION(NV,NV)          :: JJ, STM1
                                            
    integer                  :: i, j, J1, J2, J3, k, l, m, n, o
    !                                1234567890123456789012345678901234567890
    CHARACTER(LEN=40)        :: PATH='MODULE_INTEGRATION,RUNGE_KUTTA_DRIVER:'
    CHARACTER(LEN=20)        :: FMT1='(/A/,i5,3x,f12.5)'
                    
!======================================================================
!___ALLOCATE INTEGRATION VARIABLES
!======================================================================    
  
!======================================================================
!.....INITIALIZE NTUPLE SCALAR SUB-COMPONENTS

        time = t0                           !INITIALIZE INTEGRATION TIME
        step = 0                            !INITIALIZE INTEGRATION LOOP COUNTER
        done = .false.                      !ACTIVITY FLAG FOR INTEGRATION
        DELT = DT                           !STEP SIZE SET IN INIT_COND
!.....INITIALIZE INTEGRATION VARIABLES
    
        WRITE(*,*); WRITE(*,*)PATH; WRITE(*,*)"INTEGRATION HISTORY SAVED IN TIME_HISTORY(FILE=30)"
        WRITE(*,*); WRITE(*,FMT1)PATH//"STEP  TIME", STEP,TIME

        call fix_d0_dC( delt )                !INITIALIZE INTEGRATION WEIGHTS
        WRITE(ECHO_INPUT,*); WRITE(ECHO_INPUT,*) PATH; 
        WRITE(ECHO_INPUT,*)'INITIAL INTEGRATION WEIGHTS';WRITE(ECHO_INPUT,rFMT)'DELT=',DELT
        WRITE(ECHO_INPUT,rFMT)'D0=',D0; WRITE(ECHO_INPUT,rFMT)'DC=',DC
!.....EVALUATE INITIAL RK4 INTEGRATION WEIGHTS

!.....MAIN INTEGRATION LOOP
 
    DO while (.not.done)

        done = last_step( delt )          !ADJUST DELT TO HIT FINAL TIME (TF) EXACTLY

        step = step + 1                   !COUNT NUMBER OF INTEGRATION STEPS
        DO I=1,NV
            DO J=1,NV
                STM1(I,J) = X(I)%V%VPART(J)
            END DO
        END DO
        
        write(junk,*) " time = ", TIME
        write(*,*);  CALL OUT_MATRIX( junk, STM1, NV, NV, NV, 'STATE TRANSITION MATRIX = DXDX0 =')
        
        !CALL SAVE_Symplectic_Check_up_to_4TH_Order(X)  ! THIS IS TO SAVE THE SYMPLECT CHECK AT EACH TIME STEP .
        
               DO I = 1, NV
                DO J = 1, NV
                    WRITE(STM_1st_vec, *) X(I)%V%VPART(J)  ! PRINT phi 1st order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
                    DO J1 = 1, NV
                        WRITE(STM_2nd_vec, *) X(I)%T2%T2PART(J, J1)   ! PRINT phi 2nd order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
                        DO J2 = 1, NV
                            WRITE(STM_3rd_vec, *) X(I)%T3%T3PART(J, J1, J2)  ! PRINT phi 3RD order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
                            DO J3 =1, NV
                                WRITE(STM_4TH_vec, *) X(I)%T4%T4PART(J, J1,J2,J3)  ! PRINT phi 3RD order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
                            END DO  ! j3
                        END DO !J2
                    END DO !J1
                END DO  ! J
               END DO  ! I
               
        ! print in 52 time history only the state history
        write(StateTime,'(10E20.11)') (X(i)%E,i=1,NV) 
        
        CALL RK4 ( X, xC, xP, DxDT )        
        
        PRT_OR_STOP: if( (step/nprt)*nprt .eq. step .or. done ) then
             write( *,FMT1) PATH//'STEP,TIME=',step, time
            ! write(TIME_HISTORY,'(F12.5,10E20.11)') time,(G%T1(i)%E,i=1,NV)
        end if PRT_OR_STOP

        EXIT_LOOP: if( done ) then
            WRITE(*,*); WRITE(*,*) PATH; WRITE(*,*)'FINISHED IN INTEGRATION LOOP'
        else
            time = TIME + dt !FINAL INTEGRATED VALUE
        end if EXIT_LOOP

    END DO 
    
    DO I = 1, NV
        DO J = 1, NV
            WRITE(Final_Time_STM_1st_vec, *) X(I)%V%VPART(J)  ! PRINT phi 1st order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
            DO J1 = 1, NV
                WRITE(Final_Time_STM_2nd_vec, *) X(I)%T2%T2PART(J, J1)   ! PRINT phi 2nd order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
                DO J2 = 1, NV
                    WRITE(Final_Time_STM_3rd_vec, *) X(I)%T3%T3PART(J, J1, J2)  ! PRINT phi 3RD order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
                    DO J3 =1, NV
                        WRITE(Final_Time_STM_4TH_vec, *) X(I)%T4%T4PART(J, J1,J2,J3)  ! PRINT phi 3RD order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
                    END DO  ! j3
                END DO !J2
            END DO !J1
        END DO  ! J
    END DO  ! I
    
    WRITE(*,*); WRITE(*,*)PATH; WRITE(*,*)'POST INTEGRATION-CHECK FINAL TIME'
    WRITE(*,FMT1) PATH//'STEP, DELT = ', STEP, DELT
    WRITE(*,rFMT)'FINAL TIME CK: TF,TIME,TF-TIME=',TF, TIME, TF-TIME
    
   ! OUTPUT THE DATA
    WRITE(STM_HISTORY,*); WRITE(STM_HISTORY,*)PATH
    WRITE(STM_HISTORY,*)'STATE     TRANSITION MATRIX ='
    
    do i=1, NV
        do j=1,NV
            STM_1ST(I,J)= X(I)%V%VPART(J)
            do k=1,NV
                STM_2ND(I,J,K)= X(I)%T2%T2PART(J,K)
                DO L=1,NV
                    STM_3RD(I,J,K,L)= X(I)%T3%T3PART(J,K,L)
                    DO M=1,NV   
                        STM_4TH(I,J,K,L,M)= X(I)%T4%T4PART(J,K,L,M)
                    END DO
                END DO
            end do 
        end do
    end do
 
           
       ALLOCATE( DXDX0(NV, NV) )
       !DXDX0 =  X(1:NV)%V%VPART( 1:NV)
       CALL OUT_MATRIX( STM_HISTORY, STM_1ST( 1:NV, 1:NV), NV, NV, NV, 'DXDX0')
       
       if (p_order >= 2)    then
           do i=1, NV
               write(STM_HISTORY,*) 'i = ', i
               CALL OUT_MATRIX( STM_HISTORY,STM_2ND(I,1:NV,1:NV), NV, NV, NV, 'S2')
           end do
       end if
       
       if (p_order >= 3)    then
           do i=1, NV
               do j=1,NV
                   write(STM_HISTORY,*) 'i, j = ', i, j
                   CALL OUT_MATRIX( STM_HISTORY,  STM_3RD(I,J,1:NV,1:NV), NV, NV, NV, 'S3')
               end do
           end do
       end if 
       
       if (p_order >= 4)    then
           do i=1, NV
               do j=1,NV
                   do k=1,NV
                       write(STM_HISTORY,*) 'i, j, k = ', i, j, k
                       CALL OUT_MATRIX( STM_HISTORY,  STM_4TH(I,J,K,1:NV,1:NV), NV, NV, NV, 'S4')
                   end do 
               end do
           end do
       end if  
       

END SUBROUTINE RUNGE_KUTTA_DRIVER
!
!=================================================================
!
    SUBROUTINE RK4( G0, GC, GP, DGDT ) 
        IMPLICIT NONE
          
    !=======================================================================
    !
    !__PURPOSE:	Standard 4-th Order Runge-Kutta Integrator
    !           STATE IS STORED IN A GENERALIZED SCALAR FOR ALL TENSOR OPERATIONS
    !
    !__INPUTS:
    !     G0
    !     GC, GP, DGDT  PASSED TO ELIMINATE RE-ALLOCATION/DEALLOCATION INTERNALLY
    !			
    !.....THE SCALAR N-TUPLE DATA STRUCTURE IS ASSUMED TO BE
    !.....GIVEN BY:
    !                 G:= { T,  X,  P2, P3, P4 } 
    !
    !.....X     DENOTES STATE VARIABLES FOR 
    !                DXDT   = F(X,T) (OCEA Nx1)
    !
    !....S1     FIRST  ORDER STATE TRANSITION MATRIX FOR 
    !                DP2/DT = DEL(F)*P2; P2(T0) = I (DP NxN)
    !
    !....P3     SECOND ORDER STATE TRANSITION MATRIX FOR 
    !                DP3/DT = DEL(DEL(F))P2*P2+DEL(F)*P3; P3(T0) = 0 (DP NxNxN)
    !
    !....P4     THIRD ORDER STATE TRANSITION MATRIX FOR 
    !                DP4/DT = D/DT ( DEL(DEL(F))P2*P2+DEL(F)*P3) ; P4(T0) = 0 (DP NxNxNxN)
    !
    !
    !__CALLED BY: MODULE: INTEGRATION, RUNGE_KUTTA_DRIVER
    !__CALL PATH: MAINPROGRAM => RUNGE_KUTTA_DRIVER => RK4
    !
    !__NOTE: DELT NOT PASSED, DO, DC COMPUTED TO HIT FINAL TIME AND PASSED IN PROBLEM_
    !=======================================================================
    !COPYRIGHT(C) 2009 JAMES D. TURNER, ALL RIGHTS RESERVED
    !=======================================================================
    !.....ARGUMENT LIST VARIABLES
        TYPE(EB),DIMENSION(NV), INTENT(INOUT)::   G0     ! GENERALIZED SCALAR N_TUPLE 
        TYPE(EB),DIMENSION(NV), INTENT(INOUT)::   GP     ! PREDICTED   SCALAR N_TUPLE ESTIMATE
        TYPE(EB),DIMENSION(NV), INTENT(INOUT)::   GC     ! CURRENT     SCALAR N_TUPLE ESTIMATE
        TYPE(EB),DIMENSION(NV), INTENT(INOUT):: DGDT     ! CURRENT     SCALAR N_TUPLE DERIVATIVE
    !.....LOCAL VARIABLES
        REAL(DP    )              ::   TC     ! CURRENT TIME IN INTEGRATION
        INTEGER                   ::    I
    ! 
        TC = TIME 
    !
        CALL DERIVATIVE( TC, G0, DGDT )  ! k1 := function(Tc, G0)        
    !
        TC = TIME + D0(2)                ! tn+1 = tn + h/2
       
        GC = G0 + (D0(2)*DGDT)           ! Gc = G0 + (h/2)*K1
        GP = G0 + (DC(1)*DGDT)           ! Gp = G0 + (h/6)*K1

        CALL DERIVATIVE( TC, GC, DGDT )  ! k2 := function(Tc, Gc)
    !
        TC = TIME + D0(3)                ! tn+1 = tn + h/2

        GC = G0 + (D0(3)*DGDT)           ! Gc = G0 + (h/2)*K2
        GP = GP + (DC(2)*DGDT)           ! Gp = Gp + (h/3)*K2

        CALL DERIVATIVE( TC, GC, DGDT )  ! k3 := function(Tc, Gc)
    !
        TC = TIME + D0(4)                ! tn+1 = tn + h

        GC = G0 + (D0(4)*DGDT)           ! Gc = G0 + (h)*K3
        GP = GP + (DC(3)*DGDT)           ! Gp = Gp + (h/3)*K3
    !
        CALL DERIVATIVE( TC, GC, DGDT )  ! k4 := function(Tc, Gc)
    !
    !.....UPDATE TIME & STATE ESTIMATE 
    !
       ! TIME = TC
        G0   = GP + (DC(4)*DGDT)         ! G0 = Gp + (h/6)*K4 
    !
    END SUBROUTINE RK4
!
!=================================================================
!
    SUBROUTINE DERIVATIVE( TC, G, DGDT )
    !
    !__CALLED BY: MODULE: INTEGRATION, RUNGE_KUTTA_DRIVER
    !__CALL PATH: MAINPROGRAM => RUNGE_KUTTA_DRIVER => RK4 => DERIVATIVE
    !=======================================================================
    !COPYRIGHT(C) 2009 JAMES D. TURNER, ALL RIGHTS RESERVED
    !=======================================================================
    
        IMPLICIT NONE
     !.....ARGUMENT LIST VARIABLES
        REAL(DP    )     ,INTENT(IN   ):: TC
        TYPE(EB) ,DIMENSION(N_STATES)    ,INTENT(INOUT):: G          !NEED INTENT(INOUT) TO AVOID DEALLOCATION BEFORE RETURN
        TYPE(EB),DIMENSION(N_STATES)     ,INTENT(INOUT):: DGDT       !NEED INTENT(INOUT) TO AVOID DEALLOCATION BEFORE RETURN
        
    !.....LOCAL VARIABLES
        TYPE(EB    ),DIMENSION( NV ),SAVE   :: X     !OCEA DUMMY VARIABLE FOR STATE VECTOR DERIVATIVE CALCULATIONS
        TYPE(EB    ),DIMENSION( NV )        :: DXDT  !OCEA
  
        CHARACTER(LEN=31)                   :: PATH  = 'MODULE_INTEGRATION,DERIVATIVE'      
        LOGICAL(LGT),SAVE                   :: FPASS =.TRUE.
        INTEGER::I,J,K,L
       
        
    
        CALL STATE_VECTOR_DERIVATIVE ( TC, g, DgDT )
!            LOAD USER-DEFINED DIFFERENTIAL EQUATIONS AS OCEA VARIABLES                       
       
     !  call TENSOR_RATES( dxdt, G, dGdt )
!            LOAD STATE AND PARAMETER TRANSITION MATRIX DIFFERENTIAL EQUATIONS 
       
    END SUBROUTINE DERIVATIVE
!
!=================================================================
!     
   
!
!=================================================================
! 
SUBROUTINE STATE_VECTOR_DERIVATIVE ( TC, X0, DXDT )
USE PROBLEM_DATA; USE N_TUPLE; USE EB_VARIABLE_HANDLING; USE EGM2008SH200x200
IMPLICIT NONE
    !.....THIS ROUTINE LINKS THE STATE TRANSITION MATRIX CODE TO USER DERIVATIVE CALCULATOR
    !     GIVEN X0 => DXDT (OCEA), WHERE X0 = G%T2
    !..INPUTS:
    !       TC      CURRENT TIME IN THE SIMULATION
    !        G		SCALAR N-TUPLE DATA STRUCTURE
    !..OUTPUTS:
    !     DXDT      TIME DERIVATIVE OF THE SCALAR N-TUPLE DATA STRUCTURE	
    !..ASSUMPTIONS:
    !.....THE SCALAR N-TUPLE DATA STRUCTURE IS ASSUMED TO BE
    !     GIVEN BY:
    !                G:     = { T,  X,  P2, P3, Z2, Z3 } 
    !.....X     DENOTES STATE VARIABLES FOR 
    !                DXDT   = F(X,T) (OCEA Nx1)
    !....P2     FIRST  ORDER STATE TRANSITION MATRIX FOR 
    !                DP2/DT = DEL(F)*P2; P2(T0) = I (DP NxN)
    !....P3     SECOND ORDER STATE TRANSITION MATRIX FOR 
    !                DP3/DT = DEL(DEL(F))P2*P2+DEL(F)*P3; P3(T0) = 0 (DP NxNxN)
    !....Z2     FIRST  ORDER PARAMETER TRANSITION MATRIX 
    !                DZ2/DT = DEL(F)*Z2; Z2(T0) = I (DP NxN)
    !....Z3     SECOND ORDER PARAMETER TRANSITION MATRIX 
    !                DZ3/DT = DEL(DEL(F))Z2*Z2+DEL(F)*Z3; Z2(T0) = 0 (DP NxNxN)
    !=======================================================================
    !COPYRIGHT(C) 2009 JAMES D. TURNER, ALL RIGHTS RESERVED
    !=======================================================================

!.....ARGUMENT LIST VARIABLES
    REAL(DP    ),             INTENT(IN   ):: TC
    TYPE(EB    ),DIMENSION(:),INTENT(INOUT):: X0
    TYPE(EB    ),DIMENSION(:),INTENT(INOUT):: DXDT
    
    TYPE(EB    ),DIMENSION(3):: xyz
    TYPE(EB    ),DIMENSION(3):: Gxyz
    
    REAL(DP    ),DIMENSION(Max_Degree,Max_Degree):: C
    REAL(DP    ),DIMENSION(Max_Degree,Max_Degree):: S
    
!.....LOCAL VARIABLES
    real(dp    )                :: theta_m
    type(eb    )                :: rp_e, rp_e3, rp_m, rp_m3, Mue_rpe3
    real(dp    )                :: xm, ym
    TYPE(EB    )                :: P1, P2
         
    CHARACTER(LEN=44)           :: PATH = 'MODULE_INTEGRATION,STATE_VECTOR_DERIVATIVE'      
    LOGICAL(LGT)                :: FPASS = .TRUE.
    INTEGER                     :: I

!.....PROBLEM VARIABLE DEFINITIONS AND UNITS
!     ======================================
!           r_e     = Earth Radius KM
!           r_m     = Moon  Radius KM
!           v_m     = Moon inertial velocity KM/S                                 
!           omega_m = moon angular velocity RAD/S       
!           mu_e    = earth gravitational constant  KM^3/S^2
!           mu_m    = moon  gravitational constant  KM^3/S^2 
!           origin is the Earth center
!           x-axis points toward the Moon initial poistion
!           y-axis is prpenducular to the x-axis and is directed as the Moon initial inertial velocity 
! 
!.....LOAD STATE VECTOR COMPONENTS (REPLACE WITH READ STATEMENT FOR INPUT)

            IF(FPASS) THEN
               FPASS = .FALSE.
               DXDT = .CONSTANT.ZERO 
               call EGM2008_Coef( C, S)
!              INITIALIZE OCEA FIRST ORDER STATE VECTOR DERIVATIVES
            END IF
                       
            XYZ  = .CONSTANT.ZERO
            Gxyz =  .CONSTANT.ZERO
            
            rp_e  = .CONSTANT.ZERO; rp_e3 = .CONSTANT.ZERO; rp_m = .CONSTANT.ZERO
            rp_m3 = .CONSTANT.ZERO; Mue_rpe3 = .CONSTANT.ZERO
        
            rp_e    = EB_sqrt( x0(1)*x0(1) + x0(2)*x0(2) + x0(3)*x0(3) )
!                  radial distance of spacecraft from Earth        
            rp_e3    = rp_e*rp_e*rp_e
            Mue_rpe3 = mu_e/rp_e3
            dxdt(1)  = x0(4)
            dxdt(2)  = x0(5)
            dxdt(3)  = x0(6)
            
            xyz(1)   = x0(1)
            xyz(2)   = x0(2)
            xyz(3)   = x0(3)
            call EGM2008( xyz, Gxyz, C, S )
            dxdt(4)  = Gxyz(1)
            dxdt(5)  = Gxyz(2)
            dxdt(6)  = Gxyz(3)
!                 STATE DIFFERENTIAL EQUATIONS
        
        
END SUBROUTINE STATE_VECTOR_DERIVATIVE

! 
!=================================================================
!
    subroutine fix_d0_dC( delt )
    !THIS PROGRAM ADJUSTS INTEGRATION WEIGHTS
    !=======================================================================
    !COPYRIGHT(C) 2009 JAMES D. TURNER, ALL RIGHTS RESERVED
    !=======================================================================
        implicit none
        REAL(DP),INTENT(INOUT):: delt

        D0( 1 ) = 0.0D0
        D0( 2 ) = DELT/2.0e0
        D0( 3 ) = DELT/2.0e0
        D0( 4 ) = DELT

        DC( 1 ) = DELT/6.0e0
        DC( 2 ) = DELT/3.0e0
        DC( 3 ) = DELT/3.0e0
        DC( 4 ) = DELT/6.0e0

    END SUBROUTINE FIX_D0_DC
!
!=================================================================
!
    FUNCTION LAST_STEP( delt )
    ! THIS PROGRAM DETECTS THE END OF THE INTEGRATION INTERVAL AND
    ! ADJUSTS THE STEP SIZE TO EXACTLY HIT THE DESIRED FINAL TIME (TF)
    !=======================================================================
    !COPYRIGHT(C) 2009 JAMES D. TURNER, ALL RIGHTS RESERVED
    !=======================================================================
        IMPLICIT NONE
        LOGICAL                :: LAST_STEP
        REAL(DP), INTENT(INOUT):: DELT

        if( time+delt .ge. tf ) then
        
            delt      = tf - time
                !SET FINAL INTEGRATION STEPSIZE TO HIT FINAL TIME EXACTLY
                
            CALL FIX_D0_DC( DELT )
                !UPDATE THE INTEGRATION WEIGHTS
                
            last_step = .true.
        else
            last_step = .false.
        end if

    END FUNCTION LAST_STEP  

END MODULE INTEGRATION