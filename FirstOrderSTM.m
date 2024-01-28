function dX=FirstOrderSTM(~,Varin,mu) 


dR = zeros(6,1);
x = Varin(1:6);
VecPhi = Varin(7:42);
r = sqrt(x(1)^2+x(2)^2+x(3)^2);
r3 = r^3;

% State Vector (1 through 6)
dR(1) = x(4) ;
dR(2) = x(5) ;
dR(3) = x(6) ;
dR(4) = -(mu/r3)*x(1);
dR(5) = -(mu/r3)*x(2);
dR(6) = -(mu/r3)*x(3);

% First order STM (7 through 42)
dFDX41 = (3*mu*x(1)*x(1)/r^5) - mu/r3;
dFDX42 =  3*mu*x(1)*x(2)/r^5;
dFDX43 =  3*mu*x(1)*x(3)/r^5;

dFDX51 =  3*mu*x(1)*x(2)/r^5;
dFDX52 = (3*mu*x(2)*x(2)/r^5) - mu/r3;
dFDX53 =  3*mu*x(2)*x(3)/r^5;

dFDX61 =  3*mu*x(1)*x(3)/r^5;
dFDX62 =  3*mu*x(2)*x(3)/r^5;
dFDX63 = (3*mu*x(3)*x(3)/r^5) - mu/r3;

% DFDX
dFdX   = [zeros(3), eye(3);...
      dFDX41 dFDX42 dFDX43 zeros(1,3);...
      dFDX51 dFDX52 dFDX53 zeros(1,3);...
      dFDX61 dFDX62 dFDX63 zeros(1,3)];

% Calculate DX with dR and Phi*dFdX (1 through 42)
phi = reshape(VecPhi,6,6);             %Reshape from 36X1 Vector to 6X6

dphi = dFdX*phi;                       % Delta Phi at each dFdX 
dX = [dR(1:6);reshape(dphi,36,1)];     %[dR;First order STM]



end