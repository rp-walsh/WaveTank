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
Gx = @(x,z) -2*(x-pi)/c.*G(x,z);
Gz = @(x,z) -2*(z-pi)/c.*G(x,z);

% Define k and l (k = 0,+/- 1, ... l = n/2, n=0,+/- 1, ...)
k = 12; l = 12;

% Solve for discrete dispersion relation
h = dx;
b = 2-dt^2*g - dt^2*g*sin(l*h)^2/(2*cos(k*h) + 2*cos(l*h) - 4);
b_k = -2*h*dt^2*g*sin(l*h)^2*sin(h*k)/(2*cos(k*h) + 2*cos(l*h) - 4)^2;
b_l = -dt^2*g*(2*h*sin(l*h)*cos(l*h)/(2*cos(k*h) + 2*cos(l*h) - 4) + ...
               2*h*sin(l*h)^3/(2*cos(k*h) + 2*cos(l*h) - 4)^2);
r = roots([1,-b,1]);
tmp = max(real(log(r)/(1i*dt)));
omega = @(k,l) tmp;
w_k = 1/(2i*dt*exp(1i*tmp*dt) - 1i*dt*b)*b_k;
w_l = 1/(2i*dt*exp(1i*tmp*dt) - 1i*dt*b)*b_l;
Cx = w_k;
Cz = w_l;

% Define Fourier type initial data as a function of wavenumbers k and l
sTmp = @(k,l,x,z,t) (1i*k^2*G(x,z) - 2*k^2/omega(k,l)*Cx*Gx(x,z)...
                     - (2*k^2/omega(k,l)*Cz + k^2/l)*Gz(x,z))...
                     .*exp(1i*(k*x + l*z - omega(k,l)*t));
%sTmp = @(k,l,x,z,t) (1i*k^2*G(x,z) - (2*k^2*l^2*omega(k,l)^3)/(omega(k,l)*k^3*g)*Gx(x,z)...
%                     + ((2*k^2*l*omega(k,l)^3)/(omega(k,l)*k^2*g) - k^2/l)*Gz(x,z))...
%                     .*exp(1i*(k*x + l*z - omega(k,l)*t)); 
wTmp = @(k,l,x,z,t) (-k^2*omega(k,l)*G(x,z) - (1i*k^2)*Cx*Gx(x,z)...
                     - ((1i*k^2)*Cz + 1i*k^2*omega(k,l)/l)*Gz(x,z))...
                     .*exp(1i*(k*x + l*z - omega(k,l)*t));       
%wTmp = @(k,l,x,z,t) (-k^2*omega(k,l)*G(x,z) - (1i*k^2*l^2*omega(k,l)^3)/(k^3*g)*Gx(x,z)...
%                     + ((1i*k^2*l*omega(k,l)^3)/(k^2*g) - 1i*k^2*omega(k,l)/l)*Gz(x,z))...
%                     .*exp(1i*(k*x + l*z - omega(k,l)*t));       
uTmp = @(k,l,x,z,t) (k*l*omega(k,l)*G(x,z) + ...
                     (1i*k*l*Cx + 1i*l*omega(k,l))*Gx(x,z) - ...
                     1i*k*l*Cz*Gz(x,z))...
                     .*exp(1i*(k*x + l*z - omega(k,l)*t));       
%uTmp = @(k,l,x,z,t) (k*l*omega(k,l)*G(x,z) + ((1i*k*l^3*omega(k,l)^3)/(k^3*g) ...
%                                              + 1i*l*omega(k,l))*Gx(x,z)...
%                     - (1i*k*l^2*omega(k,l)^3)/(k^2*g)*Gz(x,z))...
%                     .*exp(1i*(k*x + l*z - omega(k,l)*t));       


% Define initial data for u,w,rho and BC's for P.
% Remove l or -l solution for one-way wave solution
sInit = @(x,z) 2*real(sTmp(k,l,x,z,0) + 0.*sTmp(k,-l,x,z,0));
wInit = @(x,z) 2*real(wTmp(k,l,x,z,0) + 0.*wTmp(k,-l,x,z,0));
uInit = @(x,z) 2*real(uTmp(k,l,x,z,0) + 0.*uTmp(k,-l,x,z,0));
rho = @(x,z,s) -s;
PBCT = @(x,t) zeros(size(x)); 
PBCB = @(x,t) zeros(size(x)); 
