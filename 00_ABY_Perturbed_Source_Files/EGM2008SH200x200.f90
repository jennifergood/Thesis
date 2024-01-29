MODULE EGM2008SH200x200
USE PROBLEM_DATA
USE N_TUPLE
USE EB_VARIABLE_HANDLING
USE OCEA_IO
USE OCEA_MATRIX_LIB
USE OCEA_MATHLIB
IMPLICIT NONE

  
CONTAINS
!
!=================================================================
!
SUBROUTINE EGM2008_Coef( C, S)
USE PROBLEM_DATA
    IMPLICIT NONE
    ! LOAD the Spherical Harmonic Coeff 
	! Compute gravity in ECEF coordinates
    REAL(DP    ),INTENT(OUT   ),DIMENSION(Max_Degree,Max_Degree):: C
    REAL(DP    ),INTENT(OUT   ),DIMENSION(Max_Degree,Max_Degree):: S

     read(18,*) C
     read(19,*) S

END SUBROUTINE EGM2008_Coef

SUBROUTINE EGM2008( XYZ, Gxyz, C, S )
USE PROBLEM_DATA
USE N_TUPLE
USE EB_VARIABLE_HANDLING
USE OCEA_IO
USE OCEA_MATRIX_LIB
! USE OCEA_MATHLIB 
    IMPLICIT NONE
	! determine radius of this thread's node
	TYPE(EB) ::  r 
	TYPE(EB) ::  phic  
	TYPE(EB) ::  lambda    
	TYPE(EB) ::  slambda       
	TYPE(EB) ::  clambda      
	TYPE(EB), DIMENSION(DEG+1 ) ::  smlambda
	TYPE(EB), DIMENSION(DEG+1 ) ::  cmlambda
    ! TYPE(EB) ::  XYZ 
    
    TYPE(EB ),INTENT(IN    ),DIMENSION(: ):: XYZ
    TYPE(EB ),INTENT(OUT   ),DIMENSION(: ):: Gxyz
    
    REAL(DP    ),INTENT(INOUT   ),DIMENSION(Max_Degree,Max_Degree):: C
    REAL(DP    ),INTENT(INOUT   ),DIMENSION(Max_Degree,Max_Degree):: S
    
    TYPE(EB), DIMENSION(DEG+3, deg+3 ) ::  pol
    real(dp), DIMENSION(DEG+3, deg+3 ) ::  scaleFactor
    INTEGER I, M, J, K 
    
    r       = .CONSTANT.ZERO
    phic    = .CONSTANT.ZERO
    lambda  = .CONSTANT.ZERO
    slambda = .CONSTANT.ZERO
    clambda = .CONSTANT.ZERO
    
    !R%E       = 0.0D0
    !phic%E    = 0.0D0
    !lambda%E  = 0.0D0
    !slambda%E = 0.0D0
    !clambda%E = 0.0D0


    !r%V = 1.0D0
    !phic%V = 1.0D0
    !lambda%V = 1.0D0
    !slambda%V = 1.0D0
    !clambda%V = 1.0D0

    smlambda = .CONSTANT.ZERO
    Cmlambda = .CONSTANT.ZERO
    !DO I = 1, DEG+1
    !    smlambda(I)%V%VPART(I) = ONE
    !    Cmlambda(I)%V%VPART(I) = ONE
    !END DO
    
    POL  = .CONSTANT.ZERO
    Gxyz = .CONSTANT.ZERO
    !
   ! DO I = 1, DEG+3
   !    DO J = 1, DEG+3  
   !       DO K =1, DEG+1  
   !     POL(I,J)%V%VPART(K) = ONE
   !     END DO
   ! END DO
   ! END DO
    

	! Compute geocentric radius
    R = EB_sqrt( xYZ(1)*xYZ(1) + xYZ(2)*xYZ(2) + xYZ(3)*xYZ(3) )  
    ! Compute geocentric latitude
    phic  = EB_asin( xYZ(3) / r )
    !Compute lambda                                                           
    lambda  = EB_EB_atan2( xYZ(2), xYZ(1) )

    slambda = EB_sin(lambda)
    clambda = EB_cos(lambda)
          
    smlambda(1)%E = 0.0d0
    cmlambda(1)%E = 1.0d0
    smlambda(2)   = slambda
    cmlambda(2)   = clambda
    !smlambda(2)%E = slambda%E
    !cmlambda(2)%E = clambda%E

    DO M = 3, DEG+1
        smlambda(m) = 2.0D0*clambda*smlambda(m-1) - smlambda(m-2)
        cmlambda(m) = 2.0D0*clambda*cmlambda(m-1) - cmlambda(m-2)
    END DO
    
    ! Compute normalized associated legendre polynomials
    CALL loc_gravLegendre( phic, scaleFactor, Pol )
    call loc_gravityPCPF( XYZ, Pol, C, S, smlambda, cmlambda, r, scaleFactor, Gxyz )
END SUBROUTINE EGM2008

SUBROUTINE loc_gravLegendre( phi, scaleFactor, POL )
USE PROBLEM_DATA
USE N_TUPLE
USE EB_VARIABLE_HANDLING
USE OCEA_IO
USE OCEA_MATRIX_LIB 
    IMPLICIT NONE
    ! loc_GRAVLEGENDRE internal function computing normalized associated 
    ! legendre polynomials, P, via recursion relations for spherical harmonic gravity 
	INTEGER k, p, MM, NN
    real(dp) :: m, n
    TYPE(EB) ::  cphi  
	TYPE(EB) ::  sphi 
    
    TYPE(EB),INTENT(INOUT   )::  phi
    
    TYPE(EB),INTENT(INOUT   ),DIMENSION(DEG+3, deg+3 ) ::  pol
    real(dp),INTENT(INOUT   ),DIMENSION(DEG+3, deg+3 ) ::  scaleFactor
    


    cphi = EB_cos(0.5d0*Pi - phi);
    sphi = EB_sin(0.5d0*Pi - phi);
	! Seeds for recursion formula
    POL(1,1)%e = 1.0d0            ! n = 0, m = 0;
    scaleFactor(1,1) = 0.0d0
    POL(2,1) = sqrt(3.0d0)*cphi ! n = 1, m = 0;
    scaleFactor(2,1)  = 1.0d0
    POL(2,2) = sqrt(3.0d0)*sphi ! n = 1, m = 1;
    scaleFactor(2,2) = 0.0d0
    
    DO NN = 2, DEG+2
		K = NN + 1
        n = nn
        DO Mm = 0, N
            m = mm
            p = mm + 1;
            ! Compute normalized associated legendre polynomials, P, via recursion relations 
            ! Scale Factor needed for normalization of dUdphi partial derivative
            if (n == m) THEN 
                POL(k,k) = (sqrt(2.0d0*n+1.0d0)/sqrt(2.0d0*n))*sphi*POL(k-1,k-1)
                scaleFactor(k,k) = 0.0d0
			!END IF
            else if (m == 0) THEN
                POL(k,p) = (sqrt(2.0d0*n+1.0d0)/n)*(sqrt(2.0d0*n-1.0d0)*cphi*POL(k-1,p) - ((n-1.0d0)/sqrt(2.0d0*n-3.0d0))* POL(k-2,p) )
                scaleFactor(k,p) = sqrt( (n+1.0d0)*(n)/2.0d0)
            !END IF
            else  
                POL(k,p) = (sqrt(2.0d0*n+1.0d0)/(sqrt(n+m)*sqrt(n-m)))*(sqrt(2.0d0*n-1.0d0)*cphi*POL(k-1,p) - sqrt(n+m-1.0d0)*(sqrt(n-m-1.0d0)/sqrt(2.0d0*n-3.0d0))*POL(k-2,p) )
                scaleFactor(k,p) = sqrt( (n+m+1.0d0)*(n-m));
            END IF 
        END DO !M
    END DO ! n
END SUBROUTINE loc_gravLegendre

    
SUBROUTINE loc_gravityPCPF( p, POL, C, S, smlambda, cmlambda, r, scaleFactor, Gxyz )
  USE PROBLEM_DATA
  USE N_TUPLE
  USE EB_VARIABLE_HANDLING
  USE OCEA_IO
  USE OCEA_MATRIX_LIB
  USE OCEA_MATHLIB
    IMPLICIT NONE
	! loc_GRAVITYPCPF internal function computing gravity in planet-centered
    ! planet-fixed (PCEF) coordinates using PCPF position, desired
    ! degree/order, normalized associated legendre polynomials, normalized
    ! spherical harmonic coefficients, trigonometric functions of geocentric
    ! latitude and longitude, planetary constants, and radius to center of
    ! planet. Units are MKS.

	INTEGER k, j, I
	REAL(DP)  :: K_fac, M_fac
    TYPE(EB) ::  radu, rRatio, rRatio_n, dUdrSumN, dUdphiSumN, dUdlambdaSumN, dUdrSumM, dUdphiSumM, dUdlambda, dUdr,  dUdphi,  dUdlambd, dUdlambdaSumM
    TYPE(EB) ::  X, Y, Z 
    
    TYPE(EB) :: DUMMY1, DUMMY2, DUMMY3

        
    TYPE(EB),INTENT(INOUT   ),DIMENSION(DEG+3, deg+3 ) ::  pol
    real(dp),INTENT(INOUT   ),DIMENSION(DEG+3, deg+3 ) ::  scaleFactor
    
    TYPE(EB ),INTENT(IN   ),DIMENSION(: ):: P
    TYPE(EB ),INTENT(OUT   ),DIMENSION(: ):: Gxyz
    
    REAL(DP    ),INTENT(INOUT   ),DIMENSION(Max_Degree,Max_Degree):: C
    REAL(DP    ),INTENT(INOUT   ),DIMENSION(Max_Degree,Max_Degree):: S
    
    
	TYPE(EB),INTENT(IN   ), DIMENSION(: ) ::  cmlambda
    TYPE(EB),INTENT(IN   ), DIMENSION(: ) ::  smlambda
    TYPE(EB),INTENT(IN   ) ::  r 
    INTEGER :: N, M 


    dUdrSumN        = .CONSTANT.ZERO
	dUdphiSumN      = .CONSTANT.ZERO
	dUdlambdaSumN   = .CONSTANT.ZERO
	dUdrSumM        = .CONSTANT.ZERO
	dUdphiSumM      = .CONSTANT.ZERO
	dUdlambdaSumM   = .CONSTANT.ZERO
	dUdr            = .CONSTANT.ZERO
	dUdphi          = .CONSTANT.ZERO
	dUdlambda       = .CONSTANT.ZERO
	
	! initialize summation of gravity in radial coordinates
	dUdrSumN%e      = 1.0d0
	dUdphiSumN%e    = 0.0d0
	dUdlambdaSumN%e = 0.0d0

	dUdrSumM%e      = 0.0d0
	dUdphiSumM%e    = 0.0d0
	dUdlambdaSumM%e = 0.0d0

	dUdr%e      = 0.0d0
	dUdphi%e    = 0.0d0
	dUdlambda%e = 0.0d0
 
    
     !XYZ = .CONSTANT.ZERO
     !DO I = 1, 3
     !    xyz(I)%V%VPART(I) = ONE
     !END DO
     
    x = .CONSTANT.ZERO; y = .CONSTANT.ZERO; z = .CONSTANT.ZERO
   ! DO I = 1, NV
   !      X%V%VPART(I) = ONE
   !      Y%V%VPART(I) = ONE
   !      Z%V%VPART(I) = ONE
   ! END DO
            
    x = p(1)
    y = p(2)
    z = p(3)
    
    radu%e= 0.0d0
    radu = r;
	rRatio = REQ/radu;
	rRatio_n = rRatio;
	! summation of gravity in radial coordinates
    DO N = 2, DEG
        k = n+1;
        rRatio_n = rRatio_n*rRatio
        dUdrSumM      = .CONSTANT.ZERO
        dUdphiSumM    = .CONSTANT.ZERO
        dUdlambdaSumM = .CONSTANT.ZERO
        !dUdrSumM%E      = 0.0D0
        !dUdphiSumM%E    = 0.0D0
        !dUdlambdaSumM%E = 0.0D0
        DO M = 0, N
            j = m+1;
            M_fac   = M;                 ! DEFINES TO CHANGE M FROM INTEGER TO DP
			dUdrSumM      = dUdrSumM + POL(k,j) *(C(k,j)*cmlambda(j) + S(k,j)*smlambda(j)); 
            dUdphiSumM    = dUdphiSumM + ((POL(k,j+1)*scaleFactor(k,j)) - z/(EB_SQRT(x*x + y*y))*M_fac*POL(k,j))*(C(k,j)*cmlambda(j) + S(k,j)*smlambda(j)); 
            dUdlambdaSumM = dUdlambdaSumM + M_fac*POL(k,j)*(S(k,j)*cmlambda(j) - C(k,j)*smlambda(j));
        END DO 
        K_fac         = k;                          ! DEFINES TO CHANGE K FROM INTEGER TO DP
        dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n*K_fac;
        dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
        dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;
    END DO
    ! gravity in spherical coordinates
    dUdr      = -MU/(radu*radu)*dUdrSumN ;
	dUdphi    =  MU/radu*dUdphiSumN ;
	dUdlambda =  MU/radu*dUdlambdaSumN ;
		
	! gravity in ECEF coordinates

    Gxyz(1) = ((1.0d0/radu)*dUdr - (z/(radu*radu*EB_sqrt(x*x + y*y)))*dUdphi)*x - (dUdlambda/(x*x + y*y))*y; 
    Gxyz(2) = ((1.0d0/radu)*dUdr - (z/(radu*radu*EB_sqrt(x*x + y*y)))*dUdphi)*y + (dUdlambda/(x*x + y*y))*x; 
    Gxyz(3) = (1.0d0/radu)*dUdr*z + ((EB_sqrt(x*x + y*y))/(radu*radu))*dUdphi;
	
END SUBROUTINE loc_gravityPCPF
      
END MODULE EGM2008SH200x200