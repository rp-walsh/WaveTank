%% Main driving code for WaveTank project
clear all; close all;

%% Define input file
%input = @input_FourierSoln;
%input = @input_gauss;
input = @input_WavePacket;
%input = @input_JumpRho;

%% Define spatial mesh parameters
Lx = 2*pi;
Lz = 2*pi;
m = 2^(6);
n = 2^(6)+2;

%% Create temporal computing mesh
tFinal = 10.0;
dt = tFinal/2^(6);

%% Turn plotting on/off
vis = true;

%% Call solver
[x,z,dx,dz,time,u,w,P,s,rho] = ToySolver(input,Lx,Lz,m,n,tFinal,dt,vis);

disp(['Inf norm of u: ' num2str(max(u(:)))])