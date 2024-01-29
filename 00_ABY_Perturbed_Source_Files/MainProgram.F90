PROGRAM MAIN
USE PROBLEM_DATA
USE OCEA_IO
USE EB_VARIABLE_HANDLING
USE N_TUPLE
USE INTEGRATION
IMPLICIT NONE
!.....THIS PROGRAM USES OCEA-BASED COMPUTATIONAL DIFFERENTIATION FOR GENERATING
!.....STATE AND/OR PARAMETER STATE TRANSITION TENSORS.  THE SENSITIVITY MODEL
!.....USES A GENERALIZED SCALAR N-TUPLE VARIABLE FOR PROCESSING THE STATE AND
!.....SENSITIVITY MODELS.  THE SCALAR N-TUPLE DATA STRUCTURE IS ASSUMED TO BE
!.....GIVEN BY:
!
!           G      = { T,  X,  S1, S2, S3 } 
!
!.....X     DENOTES STATE VARIABLES             : &
!           DXDT   = F(X,T) (OCEA)
!
!....S1     FIRST  ORDER STATE TRANSITION MATRIX: &
!           DS1/DT = DEL(F)*S1; S1(T0) = I
!
!....S2     SECOND ORDER STATE TRANSITION MATRIX: &
!           DS2/DT = DEL(DEL(F))*S1*S1+DEL(F)*S2; S2(T0) = 0
!
!....S3     THIRD  ORDER STATE TRANSITION MATRIX; S3(T0) = 0
!           DS3/DT = DEL(DEL(DEL(F)))*S1*S1*S1 + 2*DEL(DEL(F))*S2*S1 + &
!                    DEL(DEL(F)*S1*S2 + DEL(F)*S3

!.....DESIGN NOTE:
!           PARENTHESES REQUIRED FOR LINES CONTAINING TW0 OR MORE OPERATORS
!           I.E. REPLACE  Y = X + A*B  WITH Y = X + (A*B)
!
!           NEED ALLOCATED VARIABLES DECLARED AS INTENT(INOUT) TO BY-PASS DEALLOCATION
!
!           DXDT = F DEFINED IN MODULE:INTEGRATION, SUB.R: DERIVATIVE (OCEA)
!==========================================================================
!COPYRIGHT(C) 2009  JAMES D. TURNER ALL RIGHTS RESERVED
!AMMENDED: 2011
!==========================================================================
!=======================================================
TYPE(EB    ),DIMENSION(N_STATES):: X0
TYPE(EB    ),DIMENSION(N_STATES):: X
TYPE(TENSOR)                    :: G

!     CALL OPEN_IO
     open (18, file='C_200by200.txt')
     open (19, file='S_200by200.txt')
     
! 
     CALL INIT_COND( X0 )  
     
     X(:)  = X0(:)
                
     CALL RUNGE_KUTTA_DRIVER( X )
     
 !    CALL CLOSE_IO
     
     WRITE(*,*)'FINISHED IN MAIN'


END PROGRAM MAIN