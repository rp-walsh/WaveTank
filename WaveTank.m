%% Main driving code for WaveTank project

%% Define input file
input = @input_PerXNeuY;

%% Define spatial mesh parameters
Lx = 2*pi;
Lz = 2*pi;
m = 2^(6);
n = 2^(6)+2;

%% Create temporal computing mesh
tFinal = 10.0;
dt = tFinal/2^(5);

%% Turn plotting on/off
vis = true;

%% Call solver
[x,z,dx,dz,u,w,P,rho] = ToySolver(input,Lx,Lz,m,n,tFinal,dt,vis);