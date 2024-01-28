function dR = SecondOrderSTM(~,var,mu)
x = var(1:6);
dR = zeros(6+6*6+6*6*6,1);
x1 = x(1); x2 = x(2); x3 = x(3);
r = sqrt(x(1)^2+x(2)^2+x(3)^2);
VecPhi  = var(7:42);
VecPhi2 = var(43:258);
n = length(x);
r3 = r^3;

%% STATE VECTOR
dR(1) = x(4) ;
dR(2) = x(5) ;
dR(3) = x(6) ;
dR(4) = -(mu/r3)*x1;
dR(5) = -(mu/r3)*x2;
dR(6) = -(mu/r3)*x3;

%% 1st STM & dF/dx
dFDX41 = (3*mu*x1*x1/r^5) - mu/r3;
dFDX42 =  3*mu*x1*x2/r^5;
dFDX43 =  3*mu*x1*x3/r^5;

dFDX51 =  3*mu*x1*x2/r^5;
dFDX52 = (3*mu*x2*x2/r^5) - mu/r3;
dFDX53 =  3*mu*x2*x3/r^5;

dFDX61 =  3*mu*x1*x3/r^5;
dFDX62 =  3*mu*x2*x3/r^5;
dFDX63 = (3*mu*x3*x3/r^5) - mu/r3;

dFdX   = [zeros(3), eye(3);...
          dFDX41 dFDX42 dFDX43 zeros(1,3);...
          dFDX51 dFDX52 dFDX53 zeros(1,3);...
          dFDX61 dFDX62 dFDX63 zeros(1,3)];
      
Phi      = reshape(VecPhi,6,6);
Phidot   = dFdX*Phi;
dR(7:42) = Vec(Phidot);

%% Calculate 2nd order 

ddfdxx        = zeros(6,6,6);
ddfdxx(4,1,1) = 9*mu*x1/r^5 - 15*mu*x1*x1*x1/r^7;
ddfdxx(5,1,1) = 3*mu*x2/r^5 - 15*mu*x1*x1*x2/r^7;
ddfdxx(6,1,1) = 3*mu*x3/r^5 - 15*mu*x1*x1*x3/r^7;
ddfdxx(4,2,1) = 3*mu*x2/r^5 - 15*mu*x1*x1*x2/r^7;
ddfdxx(5,2,1) = 3*mu*x1/r^5 - 15*mu*x1*x2*x2/r^7;
ddfdxx(6,2,1) =- 15*mu*x1*x2*x3/r^7;
ddfdxx(4,3,1) = 3*mu*x3/r^5 - 15*mu*x1*x1*x3/r^7;
ddfdxx(5,3,1) =- 15*mu*x1*x2*x3/r^7;
ddfdxx(6,3,1) = 3*mu*x1/r^5 - 15*mu*x1*x3*x3/r^7;
ddfdxx(4,1,2) = 3*mu*x2/r^5 - 15*mu*x1*x1*x2/r^7;
ddfdxx(5,1,2) = 3*mu*x1/r^5 - 15*mu*x1*x2*x2/r^7;
ddfdxx(6,1,2) =- 15*mu*x1*x2*x3/r^7;
ddfdxx(4,2,2) = 3*mu*x1/r^5 - 15*mu*x1*x2*x2/r^7;
ddfdxx(5,2,2) = 9*mu*x2/r^5 - 15*mu*x2*x2*x2/r^7;
ddfdxx(6,2,2) = 3*mu*x3/r^5 - 15*mu*x2*x2*x3/r^7;
ddfdxx(4,3,2) =- 15*mu*x1*x2*x3/r^7;
ddfdxx(5,3,2) = 3*mu*x3/r^5 - 15*mu*x2*x2*x3/r^7;
ddfdxx(6,3,2) = 3*mu*x2/r^5 - 15*mu*x2*x3*x3/r^7;
ddfdxx(4,1,3) = 3*mu*x3/r^5 - 15*mu*x1*x1*x3/r^7;
ddfdxx(5,1,3) =- 15*mu*x1*x2*x3/r^7;
ddfdxx(6,1,3) = 3*mu*x1/r^5 - 15*mu*x1*x3*x3/r^7;
ddfdxx(4,2,3) =- 15*mu*x1*x2*x3/r^7;
ddfdxx(5,2,3) = 3*mu*x3/r^5 - 15*mu*x2*x2*x3/r^7;
ddfdxx(6,2,3) = 3*mu*x2/r^5 - 15*mu*x2*x3*x3/r^7;
ddfdxx(4,3,3) = 3*mu*x1/r^5 - 15*mu*x1*x3*x3/r^7;
ddfdxx(5,3,3) = 3*mu*x2/r^5 - 15*mu*x2*x3*x3/r^7;
ddfdxx(6,3,3) = 9*mu*x3/r^5 - 15*mu*x3*x3*x3/r^7;

for i=1:n                                           % Calculate Phi 2nd order
    for j=1:n
        for k=1:n
            Phi2(i,j,k) = VecPhi2(k+6*(j-1)+36*(i-1));
        end
    end
end
clear i;clear j;clear k;


for j=1:n                                           % Phi dot = D2DX2 + DFDX
    for k=1:n
        T1_sum = zeros(n,1); T2_sum  = zeros(n,1);
        for j1 =1:n
            for j2=1:n
                T1_sum = T1_sum + ddfdxx(:,j1,j2)*Phi(j1,j)*Phi(j2,k);
            end
            T2_sum = T2_sum + dFdX(:,j1)*Phi2(j1,j,k);
        end
        Phidot2(:,j,k) = T1_sum + T2_sum ;
    end
end



clear i;clear j;clear k;                            % Phi 2nd order into dR Tensor 
for i=1:n
    for j=1:n
        for k=1:n
            dR(k+6*(j-1)+36*(i-1)+42) = Phidot2(i,j,k);
        end
    end
end

for i=1:length(x)
    Phi2(:,:,i) = reshape(VecPhi2(36*(i-1)+1:36*i),6,6);    % Reshape Vector (43-258) into Phi2 Tensor matrix 
end


for ii=1:length(x)
    Phidot2(:,:,ii) = ddfdxx(:,:,ii)*Phi*Phi + dFdX*Phi2(:,:,ii);
end

for i=1:length(x)
    dR(36*(i-1)+43:36*i+42) = Vec(Phidot2(:,:,i));
end


