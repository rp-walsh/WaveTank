%% Main driving code for WaveTank project
clear all; close all;

%% Define input file
%input = @input_FourierSoln;
%input = @input_gauss;
%input = @input_WavePacket;
input = @input_AsymptoticWavePacket;
%input = @input_DiscreteDispr;
%input = @input_JumpRho;

%% Define spatial mesh parameters
Lx = 2*pi;
Lz = 2*pi;
m = 2^(8);
n = 2^(8)+2;

%% Create temporal computing mesh
tFinal = 60.0;
dt = tFinal/2^(8);

%% Turn plotting/saving on/off
vis = true;
sv = false;

%% Call solver
[x,z,dx,dz,time,u,w,P,s,rho] = ToySolver(input,Lx,Lz,m,n,tFinal,dt,vis,sv);

disp(['Inf norm of u: ' num2str(max(u(:)))])