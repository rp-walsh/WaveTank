%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This solver approximates solutions to the toy gravity wave problem
% described in README. 
% inputs:
%        input = input file for problem specificiation: 
%                 - uInit(x,z) = initial data for u
%                 - wInit(x,z) = initial data for w
%                 - sInit(x,z) = initial data for s
%                 - rho(x,z,s) = \rho as a function of position and
%                                entropy s
%                 - PBCT(x,t) = Top Neumann BC for pressure
%                 - PBCB(x,t) = Top Neumann BC for pressure
%        Lx,Lz = x and z domain boundaries respectivly
%        m,n = number of points in x and y respectivly
%        tFinal = final computing time
%        dt = time-step
%        vis = boolian visual parameter
% outputs:
%         x,z = x and z mesh coordinates (meshgrid format)
%         dx,dz = mesh spacing in x and y
%         time = final computing time (could be different from tFinal)
%         uNew,wNew = u and w solutions at t = time
%         P = P solution at t = time-dt/2
%         s_nph = s solution at t = time+dt/2
%         rho_nph = rho solution at t = time+dt/2
%
% Ray Walsh -- 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,z,dx,dz,time,uNew,wNew,P,s_nph,rho_nph] = ToySolver(input,Lx,Lz,m,n,tFinal,dt,vis)

%% Define exact solution through input file
input()

%% Create spatial computing mesh
dx = Lx/m;% corrected for periodic grid
dz = Lz/(n-2);% corrected for staggered grid
x = 0:dx:Lx-dx;
z = 0-dz/2:dz:Lz+dz/2;
[x,z] = meshgrid(x,z);
xIn = x(2:end-1,:);
zIn = z(2:end-1,:);

%% Create temporal computing mesh
N = round(tFinal/dt);

%% Initial condition for velocity
uOld = uInit(xIn,zIn);
wOld = wInit(xIn,zIn);

%% Obtain s_n+1/2 using one step of Heun's method
s0 = sInit(xIn,zIn);% define initial s
rho0 = rho(xIn,zIn,s0);% define initial rho

%% Calculate rho_z as source for initial pressure
% Assumes rho is 0 top and bottom boundary
dzrho0 = ([rho0(2:end,:); -rho0(end,:)] ...
                 -  [-rho0(1,:); rho0(1:end-1,:)] )/dz/2;

%% Define top/bottom Neumann boundary conditions for pressure
T = @(x) PBCT(x,0);
B = @(x) PBCB(x,0);

%% Solve for Pressure at t = 0
P = Poisson(-g*dzrho0,x,z,B,T);

%% Compute P_z for wTld update
P_z = (P(3:end,:) - P(1:end-2,:))/2/dz;

%% Compute \tilde{y}_{i+1} for Heun's meth. update 
wTld = wOld - (dt/2)*(P_z + g*rho0);

%% Update s using Heun's method
s_nph = s0 - (dt/4)*(wOld + wTld);

%% Clear tmp variables
clear s0 rho0 dzrho0 wTld

%% Obtain rho_n+1/2 from s_n+1/2
rho_nph = rho(xIn,zIn,s_nph);

%% Begin time-stepping
for n=1:N
    
    time = n*dt;
    
    %% Compute drho/dz at n+1/2
    
    % Assumes rho is periodic
    %dzrho_nph = ([rho_nph(2:end,:); rho_nph(1,:)] - [rho_nph(end,:); ...
    %                    rho_nph(1:end-1,:)])/dz/2;
    
    % Assumes rho is 0 top and bottom boundary
    dzrho_nph = ([rho_nph(2:end,:); -rho_nph(end,:)] ...
                 -  [-rho_nph(1,:); rho_nph(1:end-1,:)] )/dz/2;
    
    %% Solve for Pressure at n+1/2
    P = Poisson(-g*dzrho_nph,x,z,B,T);
    
    %% Compute grad P
    P_z = (P(3:end,:) - P(1:end-2,:))/2/dz;
    P_x = ([P(2:end-1,2:end) P(2:end-1,1)] - [P(2:end-1,end) P(2:end-1,1:end-1)])/2/dx;
    
    %% Step velocity forward in time
    uNew = uOld - dt*(P_x);
    wNew = wOld - dt*(P_z + g*rho_nph);
    uOld = uNew;
    wOld = wNew;
    
    %% Step entropy forward in time
    s_nph = s_nph - dt*wNew;

    %% Obtain updated density from updated entropy
    rho_nph = rho(xIn,zIn,s_nph);
    
    if vis
        %% Plot velocity and error
        figure(1)
        
        subplot(2,2,1)
        pcolor(xIn,zIn,uNew)
        title(['Computed u at time = ' num2str(time)])
        shading flat; axis equal; colorbar;

        subplot(2,2,2)
        pcolor(xIn,zIn,wNew)
        title(['Computed w at time = ' num2str(time)])
        shading flat; axis equal; colorbar;
        
        subplot(2,2,3)
        pcolor(xIn,zIn,s_nph)
        title(['Computed rho at time = ' num2str(time+dt/2)])
        shading flat; axis equal; colorbar;
        
        subplot(2,2,4)
        pcolor(x,z,P)
        title(['Computed P at time = ' num2str(time-dt/2)])
        shading flat; axis equal; colorbar;
        drawnow

    end
end
