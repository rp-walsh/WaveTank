%% Main driving code for WaveTank project
clear all; %close all;

%% Add subdirectories
addpath(genpath('Inputs'));
addpath(genpath('CoreCode'));

%% Define input file
%input = @input_FourierSoln;
%input = @input_gauss;
%input = @input_WavePacket;
%input = @input_AsymptoticWavePacket;
%input = @input_DiscreteDispr;
%input = @input_JumpRho;
input = @input_WtrLssCld;
%input = @input_WtrLssCldWvPckt;
%input = @input_WtrLssCldWvPcktBckgrnd2;
%input = @input_WtrLssCldBckgrnd;

%input();

%% Define spatial mesh parameters
%Lx = 2*pi;
%Lz = 4*pi;
Lx = 1;
Lz = 1;
m = 2^(7);
n = 2^(7)+2;

%% Create temporal computing mesh
%tFinal = 50.0;
%dt = tFinal/2^(8);
%dt = tFinal/2^(7);
tFinal = 35.0;
dt = tFinal/2^(9);

%% Turn plotting/saving on/off
vis = true;
sv = false;
err = false;

%% Call solver
%[x,z,dx,dz,time,u,w,P,s,rho] = ToySolverBckgrnd(input,Lx,Lz,m,n,tFinal,dt,vis,sv,err);
[x,z,dx,dz,time,u,w,P,s,rho] = ToySolver(input,Lx,Lz,m,n,tFinal,dt,vis,sv,err);

disp(['Inf norm of u: ' num2str(max(u(:)))])