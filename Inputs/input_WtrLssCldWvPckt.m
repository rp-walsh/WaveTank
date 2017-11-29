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
g = 9.81; omega = @(k,l) sqrt(g*k^2/(k^2+l^2)); c = .005 ;
H = 1; % height of domain
h = 1/2; % height of the interface
cp = 2; % above interface
cm = 1; % below interface
eps = 5e-6;
bbeta = 1; % Entropy coeff 

Gauss = @(x,z) exp(-((x-.25).^2 + (z-.75).^2)/c);
Gaussx = @(x,z) -2*(x-.25)/c.*Gauss(x,z);
Gaussz = @(x,z) -2*(z-.75)/c.*Gauss(x,z);
G = @(x,z) Gauss(x,z).*(Gauss(x,z)>=1e-4);
Gx = @(x,z) Gaussx(x,z).*(Gauss(x,z)>=1e-4);
Gz = @(x,z) Gaussz(x,z).*(Gauss(x,z)>=1e-4);

% Define Fourier type initial data as a function of wavenumbers k and l
sTmp = @(k,l,x,z,t) (1i*k^2*G(x,z) - (2*k^2*l^2*omega(k,l)^3)/(omega(k,l)*k^3*g)*Gx(x,z)...
                     + ((2*k^2*l*omega(k,l)^3)/(omega(k,l)*k^2*g) - k^2/l)*Gz(x,z))...
                     .*exp(1i*(k*x + l*z - omega(k,l)*t));
wTmp = @(k,l,x,z,t) (-k^2*omega(k,l)*G(x,z) - (1i*k^2*l^2*omega(k,l)^3)/(k^3*g)*Gx(x,z)...
                     + ((1i*k^2*l*omega(k,l)^3)/(k^2*g) - 1i*k^2*omega(k,l)/l)*Gz(x,z))...
                     .*exp(1i*(k*x + l*z - omega(k,l)*t));       
uTmp = @(k,l,x,z,t) (k*l*omega(k,l)*G(x,z) + ((1i*k*l^3*omega(k,l)^3)/(k^3*g) ...
                                              + 1i*l*omega(k,l))*Gx(x,z)...
                     - (1i*k*l^2*omega(k,l)^3)/(k^2*g)*Gz(x,z))...
                     .*exp(1i*(k*x + l*z - omega(k,l)*t));       

% Define k and l (k = 0,+/- 1, ... l = n/2, n=0,+/- 1, ...)
k = 12*2*pi; l = 12*2*pi;
-sqrt(g)*k*l/(k^2+l^2)^(3/2)

% Define initial data for u,w,rho and BC's for P.
% Remove l or -l solution for one-way wave solution
sInit = @(x,z) bbeta*(z-h) + eps*2*real(sTmp(k,l,x,z,0) + 0.*sTmp(k,-l,x,z,0));
wInit = @(x,z) zeros(size(x)) + eps*2*real(wTmp(k,l,x,z,0) + 0.*wTmp(k,-l,x,z,0));
uInit = @(x,z) zeros(size(x)) + eps*2*real(uTmp(k,l,x,z,0) + 0.*uTmp(k,-l,x,z,0));
rho = @(x,z,s) -cp*s.*(s>=0) - cm*s.*(s<0);

PBkg = @(x,z,s) bbeta*g*cm*(z-h).^2./2.*(s<0) + bbeta*g*cp*(z-h).^2./2.*(s>=0);
PBCT = @(x,t) (H-h)*bbeta*g*cp*ones(size(x)); 
PBCB = @(x,t) -h*bbeta*g*cm*ones(size(x)); 
