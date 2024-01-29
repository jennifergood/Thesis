SUBROUTINE  SAVE_Symplectic_Check_up_to_4TH_Order(X)
    Use INTEGRATION
    USE problem_data
    USE OCEA_IO; USE OCEA_MATRIX_LIB
IMPLICIT NONE

    TYPE(EB), DIMENSION(NV), INTENT(INOUT)  ::   X     ! GENERALIZED SCALAR N_TUPLE                                         
    REAL(DP)                               :: SUM, sum1 , sum2 , SUM3 , SUM4                                  
    REAL(DP), DIMENSION(NV,NV)             :: Symplectic_1st_order_STM, stm1, S1_t, jj
    REAL(DP), DIMENSION(NV,NV,NV)          :: Symplectic_2nd_order_STM
    real(DP), dimension(NV,NV,NV,NV)       :: Symplectic_3rd_order_STM
    real(DP), dimension(NV,NV,NV,NV,NV)    :: Symplectic_4TH_order_STM
    integer                                :: i, j, k, l, m, n, o, CNTR
    integer                                :: j1, j2, J3
    
      ! Check 1st Order Symplectic Behavior
       JJ = 0.D0
      !  Symplectic Check for MATLAB result: Phi(t, t0)
      JJ(1,4) = 1.0D0
      JJ(2,5) = 1.0D0
      JJ(3,6) = 1.0D0
      JJ(4,1) =-1.0D0
      JJ(5,2) =-1.0D0
      JJ(6,3) =-1.0D0
      
      ! SAVE THE STM 1ST ORDER 
      DO I=1,NV
          DO J=1,NV
              STM1(I,J) = X(I)%V%VPART(J)
          END DO
      END DO
      S1_t = .T.STM1
      ! Check 1ST Order Symplectic Behavior
      Symplectic_1st_order_STM = MATMUL(MATMUL(S1_t, JJ), STM1)-JJ
      ! Check 2nd Order Symplectic Behavior
       if (p_order >= 2)    then
       DO i=1, NV
          DO l=1, NV
              DO m=1, NV
                  sum = 0.0d0
                  DO j=1, NV
                     DO k=1, NV
                        sum = sum + X(j)%T2%T2PART(i,m)*JJ(j,k)*STM1(k,l) + &
                                    STM1(j,i)*JJ(j,k)*X(k)%T2%T2PART(l,m) 
                     END DO
                  END DO
                  Symplectic_2nd_order_STM(i,l,m) = sum;
              END DO
          END DO
       END DO
       END IF
       
       
       ! Check 3RD Order Symplectic Behavior
       if (p_order >= 3)    then
       DO i=1, NV
          DO l=1, NV
              DO m=1, NV
                  DO n=1, NV
                          sum = 0.0d0
                          DO j=1, NV
                              DO k=1, NV
                                  sum = sum + X(j)%T3%T3PART(i,m,n)*JJ(j,k)*STM1(k,l) + &
                                              X(j)%T2%T2PART(i,m)*JJ(j,k)*X(k)%T2%T2PART(l,n) + &
                                              X(j)%T2%T2PART(i,n)*JJ(j,k)*X(k)%T2%T2PART(l,m) + &
                                              STM1(j,i)*JJ(j,k)*X(k)%T3%T3PART(l,m,n) 
                              END DO
                          END DO
                          Symplectic_3rd_order_STM(i,l,m,n) = sum;
                  END DO
              END DO
          END DO
       END DO
       end if  
       
       ! Check 4th Order Symplectic Behavior
       if (p_order >= 4)    then
       DO i=1, NV
          DO l=1, NV
              DO m=1, NV
                  DO n=1, NV
                      DO o=1, NV
                          sum = 0.0d0
                          DO j=1, NV
                              DO k=1, NV
                                  sum = sum + X(j)%T4%T4PART(i,m,n,o)*JJ(j,k)*STM1(k,l) + &
                                              X(j)%T3%T3PART(i,m,n)*JJ(j,k)*X(k)%T2%T2PART(l,o) + &
                                              X(j)%T3%T3PART(i,m,o)*JJ(j,k)*X(k)%T2%T2PART(l,n) + &
                                              X(j)%T2%T2PART(i,m)*JJ(j,k)*X(k)%T3%T3PART(l,n,o) + &
                                              X(j)%T3%T3PART(i,n,o)*JJ(j,k)*X(k)%T2%T2PART(l,m) + &
                                              X(j)%T2%T2PART(i,n)*JJ(j,k)*X(k)%T3%T3PART(l,m,o) + &
                                              X(j)%T2%T2PART(i,o)*JJ(j,k)*X(k)%T3%T3PART(l,m,n) + &
                                              STM1(j,i)*JJ(j,k)*X(k)%T4%T4PART(l,m,n,o) 
                              END DO
                          END DO
                          Symplectic_4th_order_STM(i,l,m,n,o) = sum;
                      END DO
                  END DO
              END DO
          END DO
       END DO
    end if  
    
    
    !=========================[PRINT OUT THE SAVED SYMPLECTIC DATA]============]
    DO I = 1, NV
        DO J = 1, NV
            WRITE(SAVE_SYMPLECTIC_STM_1ST, *) Symplectic_1st_order_STM(I,J)  ! PRINT phi 1st order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
            DO J1 = 1, NV
                WRITE(SAVE_SYMPLECTIC_STM_2ND, *) Symplectic_2nd_order_STM(I, J, J1)  ! PRINT phi 2nd order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
                DO J2 = 1, NV
                    WRITE(SAVE_SYMPLECTIC_STM_3RD, *) Symplectic_3rd_order_STM(I,J,J2,J2)  ! PRINT phi 3RD order AS TIME X VEC(PHI)' ~ 300 x (1:..:36)
                    DO J3= 1, NV
                        WRITE(SAVE_SYMPLECTIC_STM_4TH, *) Symplectic_4TH_order_STM(I,J,J1,J2,J3)  ! PRINT phi 3RD order AS TIME X VEC(PHI)' ~ 300 x (1:..:36) 
                    END DO ! J3
                END DO !J2
            END DO !J1                
        END DO  ! J
    END DO  ! I
        
       
       end SUBROUTINE