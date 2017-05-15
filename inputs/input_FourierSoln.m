%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% -Input file for ToySolver.m using periodic BCs in x and homogenous
% Neumann in y for the pressure inversion.
% -This input file contains both the problem set up as well as the exact
%  solution i.e. this input file is well suited for a convergence study.
%
% This input file MUST specify atleast the following:
%                 - uInit(x,z) = initial data for u
%                 - wInit(x,z) = initial data for w
%                 - sInit(x,z) = initial data for s
%                 - rho(x,z,s) = density as a function of position and
%                                entropy 
%                 - PBCT(x,t) = Top Neumann BC for pressure
%                 - PBCB(x,t) = Top Neumann BC for pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define gravity and dispersion relation
g = 9.81; omega = @(k,l) sqrt(g*k^2/(k^2+l^2));

% Define derived Fourier solution as a function of wavenumbers k and l
uTmp = @(k,l,x,z,t) 1i*l^2*omega(k,l)*exp(1i*(k*x + l*z + omega(k,l)*t));
wTmp = @(k,l,x,z,t) -1i*k*l*omega(k,l)*exp(1i*(k*x + l*z + omega(k,l)*t));
PTmp = @(k,l,x,z,t) 1i*k*(omega(k,l)^2-g)*exp(1i*(k*x + l*z + omega(k,l)*t));
rhoTmp = @(k,l,x,z,t) -k*l*exp(1i*(k*x + l*z + omega(k,l)*t));

% Define k and l (k and l must be even to maintain 2pi periodicity)
k = 4; l = 4;

% Define exact solution. Note that we have added l and -l solutions: this
% generates a solution that satisfies homogenuous Neumann BC's in P and
% homogenuous dirichlet for w and rho (start 0 stay 0)
uExact = @(x,z,t) 2*real(uTmp(k,l,x,z,t) + uTmp(k,-l,x,z,t));
wExact = @(x,z,t) 2*real(wTmp(k,l,x,z,t) + wTmp(k,-l,x,z,t));
PExact = @(x,z,t) 2*real(PTmp(k,l,x,z,t) + PTmp(k,-l,x,z,t));
rhoExact = @(x,z,t) 2*real(rhoTmp(k,l,x,z,t) + rhoTmp(k,-l,x,z,t));
sExact = @(x,z,t) -rhoExact(x,z,t);
PExact_z = @(x,z,t) 2*real(1i*l*PTmp(k,l,x,z,t) - 1i*l*PTmp(k,-l,x,z,t));

% Define initial data for u,w,rho and BC's for P.
uInit = @(x,z) uExact(x,z,0);
wInit = @(x,z) wExact(x,z,0);
sInit = @(x,z) sExact(x,z,0);
rho = @(x,z,s) -s;
PBCT = @(x,t) zeros(size(x)); 
PBCB = @(x,t) zeros(size(x)); 
