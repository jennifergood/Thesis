SUBROUTINE  Symplectic_Check_up_to_3rdOrder(G0)
    Use INTEGRATION
    USE problem_data
    USE OCEA_IO; USE OCEA_MATRIX_LIB
IMPLICIT NONE

    TYPE(TENSOR),INTENT(INOUT)::   G0     ! GENERALIZED SCALAR N_TUPLE                                         
    
    REAL(DP)                 :: SUM, sum1 , sum2 , SUM3 , SUM4    
                               !INTEGRATION STEP SIZE
                               
    REAL(DP), DIMENSION(NV,NV,NV,NV,NV) :: Symplectic_4th_order_STM
    
    REAL(DP), DIMENSION(NV,NV)          :: JJ, Symplectic_1st_order_STM, S1_t
    
    REAL(DP), DIMENSION(NV,NV,NV)       ::Symplectic_2nd_order_STM
    
    real(dp), dimension(nv,nv,nv,nv)    ::Symplectic_3rd_order_STM
                                            
    integer                  :: i, j, k, l, m, n, o, CNTR
    
       ! Check 1st Order Symplectic Behavior
       JJ = 0.D0
       
       !  Symplectic Check for MATLAB result: Phi(t, t0)
      JJ(1,4) = 1.0D0
      JJ(2,5) = 1.0D0
      JJ(3,6) = 1.0D0
      JJ(4,1) =-1.0D0
      JJ(5,2) =-1.0D0
      JJ(6,3) =-1.0D0
      

          
       S1_t = .T.G0%S1
       Symplectic_1st_order_STM = MATMUL(MATMUL(S1_t, JJ), G0%S1)-JJ
       CALL OUT_MATRIX( SYMPLECTIC_STM_1ST, Symplectic_1st_order_STM, NV, NV, NV)
       
       ! Check 2nd Order Symplectic Behavior
       if (p_order >= 2)    then
       DO i=1, NV
          DO l=1, NV
              DO m=1, NV
                  sum = 0.0d0
                  DO j=1, NV
                     DO k=1, NV
                        sum = sum + G0%S2(j,i,m)*JJ(j,k)*G0%S1(k,l) + &
                                    G0%S1(j,i)*JJ(j,k)*G0%S2(k,l,m)
                     END DO
                  END DO
                  Symplectic_2nd_order_STM(i,l,m) = sum;
              END DO
          END DO
       END DO
       do i=1, NV
          write(SYMPLECTIC_STM_2ND,*) 'i = ', i
          CALL OUT_MATRIX( SYMPLECTIC_STM_2ND, Symplectic_2nd_order_STM(:,:,i), NV, NV, NV)
       end do;END IF
       
       
       ! Check 3RD Order Symplectic Behavior
       if (p_order >= 3)    then
       DO i=1, NV
          DO l=1, NV
              DO m=1, NV
                  DO n=1, NV
                          sum = 0.0d0
                          DO j=1, NV
                              DO k=1, NV
                                  sum = sum + G0%S3(j,i,m,n)*JJ(j,k)*G0%S1(k,l) + &
                                              G0%S2(j,i,m)*JJ(j,k)*G0%S2(k,l,n) + &
                                              G0%S2(j,i,n)*JJ(j,k)*G0%S2(k,l,m) + &
                                              G0%S1(j,i)*JJ(j,k)*G0%S3(k,l,m,n) 
                              END DO
                          END DO
                          Symplectic_3rd_order_STM(i,l,m,n) = sum;
                  END DO
              END DO
          END DO
       END DO

       do i=1, NV
          do j=1,NV
                write(SYMPLECTIC_STM_3RD,*) 'i, j = ', i, j
                CALL OUT_MATRIX( SYMPLECTIC_STM_3RD, Symplectic_3rd_order_STM(:,:,j,i), NV, NV, NV)
          end do
       end do; end if  
       
       ! Check 4th Order Symplectic Behavior
       !if (p_order >= 4)    then
       !DO i=1, NV
       !   DO l=1, NV
       !       DO m=1, NV
       !           DO n=1, NV
       !               DO o=1, NV
       !                   sum = 0.0d0
       !                   DO j=1, NV
       !                       DO k=1, NV
       !                           sum = sum + G0%S4(j,i,m,n,o)*JJ(j,k)*G0%S1(k,l) + &
       !                                       G0%S3(j,i,m,n)*JJ(j,k)*G0%S2(k,l,o) + &
       !                                       G0%S3(j,i,m,o)*JJ(j,k)*G0%S2(k,l,n) + &
       !                                       G0%S2(j,i,m)*JJ(j,k)*G0%S3(k,l,n,o) + &
       !                                       G0%S3(j,i,n,o)*JJ(j,k)*G0%S2(k,l,m) + &
       !                                       G0%S2(j,i,n)*JJ(j,k)*G0%S3(k,l,m,o) + &
       !                                       G0%S2(j,i,o)*JJ(j,k)*G0%S3(k,l,m,n) + &
       !                                       G0%S1(j,i)*JJ(j,k)*G0%S4(k,l,m,n,o) 
       !                       END DO
       !                   END DO
       !                   Symplectic_4th_order_STM(i,l,m,n,o) = sum;
       !               END DO
       !           END DO
       !       END DO
       !   END DO
       !END DO

       !do i=1, NV
        !  do j=1,NV
        !     do k=1,NV
        !        write(STM_HISTORY,*) 'i, j, k = ', i, j, k
        !        CALL OUT_MATRIX( STM_HISTORY, Symplectic_4th_order_STM(:,:,k,j,i), NV, NV, NV, 'Symplectic_4th_order_STM')
        !     end do 
        !  end do
       !end do; end if  
       
       end SUBROUTINE