%% Input file for ToySolver2.m using periodic BCs in x and Neumann in y for
%% the pressure inversion

% Define exact solution
g = 9.81; omega = @(k,l) sqrt(g*k^2/(k^2+l^2));

uTmp = @(k,l,x,z,t) 1i*l^2*omega(k,l)*exp(1i*(k*x + l*z + omega(k,l)*t));
wTmp = @(k,l,x,z,t) -1i*k*l*omega(k,l)*exp(1i*(k*x + l*z + omega(k,l)*t));
PTmp = @(k,l,x,z,t) 1i*k*(omega(k,l)^2-g)*exp(1i*(k*x + l*z + omega(k,l)*t));
rhoTmp = @(k,l,x,z,t) -k*l*exp(1i*(k*x + l*z + omega(k,l)*t));

k = 4; l = 4;

uExact = @(x,z,t) 2*real(uTmp(k,l,x,z,t) + uTmp(k,-l,x,z,t));
wExact = @(x,z,t) 2*real(wTmp(k,l,x,z,t) + wTmp(k,-l,x,z,t));
PExact = @(x,z,t) 2*real(PTmp(k,l,x,z,t) + PTmp(k,-l,x,z,t));
rhoExact = @(x,z,t) 2*real(rhoTmp(k,l,x,z,t) + rhoTmp(k,-l,x,z,t));

PExact_z = @(x,z,t) 2*real(1i*l*PTmp(k,l,x,z,t) - 1i*l*PTmp(k,-l,x,z,t));

%PExact_x = @(x,z,t) -k^2*(omega^2-g)*exp(1i*(k*x + l*z + omega*t)) + ...
%    -k^2*(omega^2-g)*exp(-1i*(k*x + l*z + omega*t));
%PExact_z = @(x,z,t) -l*k*(omega^2-g)*exp(1i*(k*x + l*z + omega*t)) + ...
%    -l*k*(omega^2-g)*exp(-1i*(k*x + l*z + omega*t));
%rhoExact_z = @(x,z,t) -1i*k*l^2*exp(1i*(k*x + l*z + omega*t)) + ...
%    1i*k*l^2*exp(-1i*(k*x + l*z + omega*t));
