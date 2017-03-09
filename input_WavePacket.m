%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% - Input file for ToySolver.m using periodic BCs in x and homogenous
%   Neumann in y for the pressure inversion.
% - This input file contains ONLY the problem set up and NOT the exact
%   solution i.e. this input file is IS NOT well suited for a convergence
%   study. 
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
g = 9.81; omega = @(k,l) sqrt(g*k^2/(k^2+l^2)); c = 1. ;
G = @(x,z) exp(-((x-pi).^2 + (z-pi).^2)/c);

% Define Fourier type initial data as a function of wavenumbers k and l
sTmp = @(k,l,x,z,t) G(x,z).*exp(1i*(k*x + l*z + omega(k,l)*t)).*(-2*(x-pi)/c + 1i*k);
wTmp = @(k,l,x,z,t) -1i*omega(k,l)*G(x,z).*exp(1i*(k*x + l*z + omega(k,l)*t))...
       .*(-2*(x-pi)/c + 1i*k);
uTmp = @(k,l,x,z,t) 1i*omega(k,l)*G(x,z).*exp(1i*(k*x + l*z + omega(k,l)*t))...
       .*(-2*(z-pi)/c + 1i*l);

% Define k and l (k = 0,+/- 1, ... l = n/2, n=0,+/- 1, ...)
k = 6; l = 6;

% Define initial data for u,w,rho and BC's for P.
% Remove l or -l solution for one-way wave solution
sInit = @(x,z) 2*real(sTmp(k,l,x,z,0) + 0.*sTmp(k,-l,x,z,0));
wInit = @(x,z) 2*real(wTmp(k,l,x,z,0) + 0.*wTmp(k,-l,x,z,0));
uInit = @(x,z) 2*real(uTmp(k,l,x,z,0) + 0.*uTmp(k,-l,x,z,0));
rho = @(x,z,s) -s;
PBCT = @(x,t) zeros(size(x)); 
PBCB = @(x,t) zeros(size(x)); 
