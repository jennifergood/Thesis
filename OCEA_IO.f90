MODULE OCEA_IO
USE PROBLEM_DATA

!******************************************************
!  input output utilities
!******************************************************
! EB_PRINT      PRINT OCEA SCALAR VARIABLES FOR 1st-4TH ORDER
! EB_PRT        PRINT OCEA SCALAR VARIABLES FOR 2ND ORDER
! OUT_MATRIX    FORMATS MATRIX - VECTOR OUTPUTS

INTERFACE PRT_D
   MODULE PROCEDURE SCAL_PRT_D, VEC_PRT_D
END INTERFACE


CONTAINS

!=======================================================

SUBROUTINE EB_PRT( OUT, EBVar, TITLE )
    USE PROBLEM_DATA
    USE EB_VARIABLE_HANDLING
    IMPLICIT NONE
    INTEGER,          INTENT(IN):: OUT
    TYPE(EB),         INTENT(IN):: EBVar
    CHARACTER(LEN=*), INTENT(IN):: TITLE
    CHARACTER(LEN=25):: FMT
    CHARACTER(LEN= 3):: str1
    write(str1,'(i2)')nv
    FMT='(/A,/'//str1//'('//str1//'E13.5/)/)'
    WRITE(OUT,*)
    WRITE(OUT,'(A)') TITLE
    WRITE(OUT,*)
    WRITE(OUT,FMT)"SCALAR PART",EBVAR%E
    WRITE(OUT,FMT)"VECTOR PART",EBVAR%V
    WRITE(OUT,FMT)"TENSOR PART",EBVAR%T2
END SUBROUTINE EB_PRT

subroutine eb_print( out, a )
     use problem_data
     use eb_variable_handling
     implicit none
     integer, intent(in):: out
     type(eb),intent(in)::  a
     character(len=25):: fmt
     character(len= 3):: n1,n2

!.....adjust printing format statement

     write(n1,'(i2)')nv      !integer to character conversion
     write(n2,'(i2)')nv**3   !integer to character conversion
     fmt = '('//n2//'('//n1//'('//n1//'e13.5/)/))'

!.....print the buggers!!!

     write(out,*  )'function='
     WRITE(OUT,FMT)a%e

     write(out,*  )'jacobian='
     write(out,fmt)a%v

     if( deriv_order >= 2 )then
         write(out,*  )'hessian='
         write(out,fmt)a%t2
     end if

     if( deriv_order >= 3 )then
         write(out,*  )'jacobian(hessian)='
         write(out,fmt)a%t3
     end if

     if( deriv_order >= 4 )then
         write(out,*  )'hessian(hessian)='
         write(out,fmt)a%t4
     end if

end subroutine eb_print

SUBROUTINE EB_VEC_PRINT( OUT, V, TITLE )
USE EB_VARIABLE_HANDLING
IMPLICIT NONE
INTEGER              :: OUT
TYPE(EB),DIMENSION(:):: V
CHARACTER(LEN=*)     :: TITLE
character(len=2)     :: n
character(len=1)     :: lp  = '('
character(len=3)     :: n2
character(len=6)     :: fmt = "e13.5/"
character(len=8)     :: pt
character(len=50)    :: prt

write(n, '(i2)')nv
write(n2,'(i3)')nv*nv
pt  = n//fmt
prt = '(A,/'//pt//',/A,/'//n//lp//pt//'/),A,/'//n2//lp//pt//'/)/)'
WRITE(OUT,*)
WRITE(OUT,*)'========================================================='
WRITE(OUT,*)TITLE
WRITE(OUT,*)'========================================================='
WRITE(OUT,PRT)&
         'CK LINEAR EQUATION ERROR ',        v%E,&
         'CK JACOBIAN OF LINEAR EQN ERROR ', v%V,&
         'CK HESSIAN  OF LINEAR EQN ERROR ', v%T2
END SUBROUTINE EB_VEC_PRINT

!=====================================================================

!++++++++++++++JAMES D. TURNER++++++++++++++++++++++++++++++++++++++++
!
!	THIS PROGRAM PRINTS MATRIX & VECTOR OUTPUTS WITH THE FORMAT
!	STATEMENTS BEING DYNAMICALLY ADJUSTED.
!
!	INPUTS:
!		OUT:	INTEGER FOR THE OPEN OUTPUT UNIT
!		A:	MATRIX/VECTOR TO BE PRINTED (ASSUMED 6X6)
!		NA;	LEADING DIMENSION OF MATRIX A
!		
!CURRENT NUMBER OF ROWS
!		NCOL:	CURRENT NUMBER OF COLUMNS
!		TITLE:	CHARACTER STRING TO LABLE PRINTOUT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!	COPYRIGHT @ 1995, JAMES D. TURNER, ALL RIGHTS RESERVED
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
SUBROUTINE OUT_MATRIX( OUT, A, NA, NROW, NCOL)
    USE PROBLEM_DATA
    IMPLICIT NONE
!
!.....PARAMETER STATEMENTS
!
    CHARACTER(LEN= 1),PARAMETER:: C = '('
    CHARACTER(LEN= 8),PARAMETER:: D = 'E15.9/))'    !12.4
    CHARACTER(LEN=12),PARAMETER:: E = '(2A,I2,A,I2/'

!.....ARGUMENT LIST VARIABLES

    INTEGER,                    INTENT(IN):: OUT, NROW, NCOL, NA
    REAL(DP),DIMENSION(NA,NCOL),INTENT(IN):: A
 !  CHARACTER(LEN= *),          INTENT(IN):: TITLE
!
!.....LOCAL VARIABLES
!
    CHARACTER(LEN=25):: FMT
    CHARACTER(LEN= 5):: TMP
    CHARACTER(LEN= 3):: N1,N2
    INTEGER          :: I,J

    write(N1,'(i3)')nrow  !integer to character conversion
    write(N2,'(i3)')ncol  !integer to character conversion

    TMP = TRIM(ADJUSTL(N1))//C//TRIM(ADJUSTL(N2))
    FMT = E//TMP//D    !SET PRINTING FORMAT STATEMENT BASED ON INPUT
!
 WRITE(OUT,*) ((A(I, J), J=1, NCOL), I=1, NROW)
!
END SUBROUTINE OUT_MATRIX

FUNCTION SCAL_PRT_D( OUT, VAR, TITLE )
!=================================================================
!..PURPOSE:
!       THIS PROGRAM PRINTS OUT SCALAR DATA AND DEPENDENT DERIVATIVE
!       VARIABLES.
!=================================================================
!..COPYRIGHT(C) 2005 JAMES D. TURNER ALL RIGHTS RESERVED
!=================================================================
    USE EB_VARIABLE_HANDLING
    INTEGER ,        INTENT(IN):: OUT      !OUTPUT FILE FOR PRINT
    TYPE(EB),        INTENT(IN):: VAR      !VARIABLE TO BE PRINTED
    CHARACTER(LEN=*),INTENT(IN):: TITLE    !NAME OF VARIABLE
    CHARACTER(LEN=5)           :: SCAL_PRT_D
    CHARACTER(LEN=3)           :: N        !INDEX OF DEPENDENT DERIVATIVE VARIABLE
    INTEGER:: J                            !DO LOOP VARIABLES
    
    SCAL_PRT_D = 'BEGIN'
    WRITE(OUT,*)
    WRITE(OUT,'(A)')TITLE
!.....PRINT SCALAR    
    WRITE(OUT,'(A,E13.5)')TITLE//'%E =   ',VAR%E
!.....PRINT DEPENDENT DERIVATIVE VARIABLES
    DO J = 1, NV    
        IF( abs( VAR%V%VPART(J) ) >= 1.0d-6 ) &
            WRITE(OUT,'(A,I3,E13.5)')TITLE//'%V =',J,VAR%V%VPART(J)
    END DO !J
    
    SCAL_PRT_D = 'DONE '
END FUNCTION SCAL_PRT_D

FUNCTION VEC_PRT_D( OUT, VAR, TITLE )
!=================================================================
!..PURPOSE:
!       THIS PROGRAM PRINTS OUT VECTOR DATA AND DEPENDENT DERIVATIVE
!       VARIABLES.
!=================================================================
!..COPYRIGHT(C) 2005 JAMES D. TURNER ALL RIGHTS RESERVED
!=================================================================
    USE PROBLEM_DATA
    USE EB_VARIABLE_HANDLING
    INTEGER ,             INTENT(IN):: OUT      !OUTPUT FILE FOR PRINT
    TYPE(EB),DIMENSION(:),INTENT(IN):: VAR      !VARIABLE TO BE PRINTED
    CHARACTER(LEN=*)     ,INTENT(IN):: TITLE    !NAME OF VARIABLE
    CHARACTER(LEN=5)                :: VEC_PRT_D
    CHARACTER(LEN=3)                :: N        !INDEX OF DEPENDENT DERIVATIVE VARIABLE
    INTEGER::I,J                                !DO LOOP VARIABLES
    
    VEC_PRT_D = 'BEGIN'
    WRITE(OUT,*)
    WRITE(OUT,'(A)')TITLE
 !.....PRINT VECTOR ELEMENTS   
    DO I = 1, SIZE(VAR)
        WRITE(N,'(I3)')I
        WRITE(OUT,'(A,E13.5)')TITLE//'('//N//')'//'%E =   ',VAR(I)%E
 !.....PRINT DEPENDENT DERIVATIVE VARIABLES
        DO J = 1, NV    
            IF( abs( VAR(I)%V%VPART(J) ) >= 1.0d-6 ) &
                WRITE(OUT,'(A,I3,E13.5)')TITLE//'('//N//')'//'%V =',J,VAR(I)%V%VPART(J)
        END DO !J
    END DO
    
    VEC_PRT_D = 'DONE '
END FUNCTION VEC_PRT_D

END MODULE OCEA_IO

