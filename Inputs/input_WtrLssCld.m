%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% - Input file for WaveTank.m using periodic BCs in x and homogenous
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

% Define solution parameters
g = 9.81;
bbeta = 1; % Entropy coeff 
cp = 2; % above interface
cm = 1; % below interface
H = 1; % height of domain
h = 1/2; % height of the interface
%h = 0.4167; % height of the interface
%h = sqrt(2)/2; % height of the interface
%h = 0.45; % height of the interface
n = 2; % fourier number (even?)
%n = 14; % fourier number (even?)
k = n*pi; % x wave number

%% Define verticle wave strcture
alphap = @(k,omega) sqrt(g*cp*k.^2./omega.^2 - k.^2);
alpham = @(k,omega) sqrt(g*cm*k.^2./omega.^2 - k.^2);

%% Define dispersion relation
F = @(k,w) alpham(k,w).*sin(alpham(k,w)*h).*cos(alphap(k,w)*(h-H)).*...
    (g*(cp-cm)./(g*cm-w.^2)+1) - alphap(k,w).*sin(alphap(k,w)*(h-H)).*cos(alpham(k,w)*h);

omega = fzero(@(w) F(k,w),sqrt(g/2));

%% Define functions
alphap = alphap(k,omega);
alpham = alpham(k,omega);
%eps = 0.1;
%eps = .01;
%eps = .005;
eps = .05;
Ap = eps*cos(alpham*h); Am = eps*cos(alphap*(h-H));
sHatM = @(x,z,t) -alpham*Am/(cm*g-omega^2).*sin(alpham*z).*sin(k*x-omega*t);
Gmm = @(x,t) h - sHatM(x,h,t)/bbeta;

Pm = @(x,z,t) bbeta*g*cm*(z-h).^2./2 + Am*cos(alpham*z).*sin(k*x-omega*t); 
Pp = @(x,z,t) bbeta*g*cp*(z-h).^2./2 + Ap*cos(alphap*(z-H)).*sin(k*x-omega*t); 
PExact =@(x,z,t) Pm(x,z,t).*((Gmm(x,t)-z)>=0) + Pp(x,z,t).*((z-Gmm(x,t))>0);
%PBkg = @(x,z,t) bbeta*g*cm*(z-h).^2./2.*((Gmm(x,t)-z)>=0) + bbeta*g*cp*(z-h).^2./2.*((z-Gmm(x,t))>0);
PBkg = @(x,z,s) bbeta*g*cm*(z-h).^2./2.*(s<0) + bbeta*g*cp*(z-h).^2./2.*(s>=0);

Pmx = @(x,z,t) k*Am.*cos(alpham*z).*cos(k*x-omega*t);
Ppx = @(x,z,t) k*Ap.*cos(alphap*(z-H)).*cos(k*x-omega*t);
PxExact = @(x,z,t) Pmx(x,z,t).*((Gmm(x,t)-z)>=0) + Ppx(x,z,t).*((z-Gmm(x,t))>0);

Pmz = @(x,z,t) bbeta*g*cm*(z-h) - alpham*Am*sin(alpham*z).*sin(k*x-omega*t); 
Ppz = @(x,z,t) bbeta*g*cp*(z-h) - alphap*Ap*sin(alphap*(z-H)).*sin(k*x-omega*t); 
PzExact = @(x,z,t) Pmz(x,z,t).*((Gmm(x,t)-z)>=0) + Ppz(x,z,t).*((z-Gmm(x,t))>0);

Sm = @(x,z,t) bbeta*(z-h) - alpham*Am/(cm*g-omega^2).*sin(alpham*z).*sin(k*x-omega*t);
Sp = @(x,z,t) bbeta*(z-h) - alphap*Ap/(cp*g-omega^2).*sin(alphap*(z-H)).*sin(k*x-omega*t); 
%SExact =@(x,z,t) Sm(x,z,t).*heaviside(Gmm(x,t)-z) + Sp(x,z,t).*heaviside(z-Gmm(x,t));
SExact =@(x,z,t) Sm(x,z,t).*((Gmm(x,t)-z)>=0) + Sp(x,z,t).*((z-Gmm(x,t))>0);
sInit = @(x,z) SExact(x,z,0);

rhoExact =@(x,z,t) -cm*Sm(x,z,t).*heaviside(Gmm(x,t)-z) - cp*Sp(x,z,t).*heaviside(z-Gmm(x,t));
rho =@(x,z,s) -cp*s.*(s>=0) - cm*s.*(s<0);

drhodz2 = @(x,z,t) -cm*bbeta + cm*alpham.^2*Am/(cm*g-omega^2).*cos(alpham*z).*sin(k*x-omega*t);%minus
drhodz1 = @(x,z,t) -cp*bbeta + cp*alphap.^2*Ap/(cp*g-omega^2).*cos(alphap*(z-H)).*sin(k*x-omega*t);%plus
SExact =@(x,z,t) Sm(x,z,t).*((Gmm(x,t)-z)>=0) + Sp(x,z,t).*((z-Gmm(x,t))>0);

Wm = @(x,z,t) -omega/(g*cm-omega^2)*Am*alpham.*sin(alpham*z).*cos(k*x-omega*t);
Wp = @(x,z,t) -omega/(g*cp-omega^2)*Ap*alphap.*sin(alphap*(z-H)).*cos(k*x-omega*t);
wExact =@(x,z,t) Wm(x,z,t).*heaviside(Gmm(x,t)-z) + Wp(x,z,t).*heaviside(z-Gmm(x,t));
wInit = @(x,z) wExact(x,z,0);

Um = @(x,z,t) k/omega*Am*cos(alpham*z).*sin(k*x-omega*t); 
Up = @(x,z,t) k/omega*Ap*cos(alphap*(z-H)).*sin(k*x-omega*t); 
uExact =@(x,z,t) Um(x,z,t).*heaviside(Gmm(x,t)-z) + Up(x,z,t).*heaviside(z-Gmm(x,t));
uInit = @(x,z) uExact(x,z,0);

PBCT = @(x,t) (H-h)*bbeta*g*cp*ones(size(x)); 
PBCB = @(x,t) -h*bbeta*g*cm*ones(size(x)); 
