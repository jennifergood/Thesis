% -----------------------------------------------------------------------
%   File Name   : Jacobi_Integral.m
%   Compiler    : MATLAB 7.11.0 (R2010b)
%   Created by  : Ahmad Bani Younes
%   Modified by : Donghoon Kim
%   Affiliation : Aerospace Engineering Department, Texas A&M University
%   Description : It calculates Hamiltonian to validate accuract for
%                 solutions.
%   References  : Orthogonal Polynomial Approximation in Higher
%                 Higher Dimensions, p. 141, Ahmad Bani Younes 
%   Subfiles    : egm2008GPsphericalharmonic.m
%   Last Update : 2013.11.21
% -----------------------------------------------------------------------
%   < Input >
%   soln        : Set of obtained solutions
%
%   < Output >
%   H           : Hamiltonian
% -----------------------------------------------------------------------

function H = Jacobi_Integral(soln)
    
    global omega Deg
    
    rB   = soln(:,1:3);         % Position in rotating frame
    vB   = soln(:,4:6);         % Velocity in rotating frame
    term = 0.5*omega*omega.*( rB(:,1).^2 + rB(:,2).^2 );
    KE   = 0.5*( vB(:,1).^2 + vB(:,2).^2 + vB(:,3).^2 );
    V    = -egm2008GPsphericalharmonic( rB, Deg); 
    H    = V + KE - term;       % Hamiltonian    