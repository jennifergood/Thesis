function dX = EoM(~,x, mu)      % Monte Carlo
    x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); x5 = x(5); x6 = x(6);
    dX = [x4;x5;x6;
        -mu*x1/sqrt(x1^2+x2^2+x3^2)^3; 
        -mu*x2/sqrt(x1^2+x2^2+x3^2)^3; 
        -mu*x3/sqrt(x1^2+x2^2+x3^2)^3];  %ODE
end