% -----------------------------------------------------------------------
%   File Name   : FnG.m
%   Compiler    : MATLAB 7.10.0 (R2010a)
%   Created by  : Donghoon Kim
%   Affiliation : Aerospace Engineering Department, Texas A&M University
%   Description : It generates F&G solutions given time, position, and 
%                 velocity information.
%   References  : Analytical Mechanics of Space Systems, 2nd, chap. 9
%                 by Schaub and Junkins
%   Last Update : 2013.12.02
% -----------------------------------------------------------------------
%
%   < Input >
%   t0          : Current time (sec)
%   t           : Propagation time (sec)
%   r0          : Current position vector
%   v0          : Current velocity vector
%
%   < Output > 
%   r           : Propagated position vector
%   v           : Propagated velocity vector
%
% -----------------------------------------------------------------------

function [ r, v ] = FnG(t0, t, r0, v0)

    global GM

    %% Transform Row Vecto to Column Vector
    if isrow(r0); r0 = r0'; end
    if isrow(v0); v0 = v0'; end

    %% Find Mean Anomlay
    R0     = sqrt(r0'*r0);          % Magnitude of current position
    V0     = sqrt(v0'*v0);          % Magnitude of current velocity
    sigma0 = dot(r0,v0)/sqrt(GM);   % Defined 
    A      = 2/R0 - V0^2/GM;        % Reciprocal of 1/a
    a      = 1/A;                   % Semi-major axis
    M      = sqrt(GM/a^3)*(t - t0); % Mean anomaly (rad)

    %% Run Newton-Raphson Method 
    tol   = 1e-15;      % Tolerance
    itr   = 0;          % Initial iteration number
    MaxIt = 10;         % Maximum iteration number
    Ehat  = M;          % Initial guess for eccentric anomaly (rad)
    dEhat = 1;          % Initial eccentric anomaly error (rad)
    while abs(dEhat) > tol;     
        err   = M - (Ehat - (1 - R0/a)*sin(Ehat) + sigma0/sqrt(a)*(1 - cos(Ehat)));
        derr  = - 1 + (1 - R0/a)*cos(Ehat) - sigma0/sqrt(a)*sin(Ehat);
        Ehat  = Ehat - err/derr;
        dEhat = - err/derr;
        itr   = itr + 1;
        if itr > MaxIt, break, end    
    end

    %% Generate F & G Solutions
    R    = a + (R0 - a)*cos(Ehat) + sqrt(a)*sigma0*sin(Ehat);
    F    = 1 - a/R0*(1 - cos(Ehat));
    G    = (t - t0) + sqrt(a^3/GM)*(sin(Ehat) - Ehat);
    Fdot = - sqrt(GM*a)/(R*R0)*sin(Ehat);
    Gdot = 1 - a/R*(1 - cos(Ehat));
    r    = F*r0 + G*v0;
    v    = Fdot*r0 + Gdot*v0;
    
end