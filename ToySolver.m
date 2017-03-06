%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This solver approximates solutions to the toy gravity wave problem
% described in README. 
% inputs:
%        input = input file for problem specificiation: 
%                 - uInit(x,z) = initial data for u
%                 - wInit(x,z) = initial data for w
%                 - sInit(x,z,t) = initial data for s (time-dependent
%                                  to maintain conv study functionality)
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

%% Obtain s_n+1/2 (use s(0) when no exact solution)
s_nph = sInit(xIn,zIn,dt/2);

%% Obtain rho_n+1/2 from s_n+1/2
rho_nph = rho(xIn,zIn,s_nph);

%% Begin time-stepping
for n=1:N
    
    time = n*dt;
    
    %% Compute drho/dz at n+1/2
    %dzrho_nph = ([rho_nph(2:end,:); rho_nph(1,:)] - [rho_nph(end,:); ...
    %                    rho_nph(1:end-1,:)])/dz/2;% Assumes rho is periodic
    dzrho_nph = ([rho_nph(2:end,:); -rho_nph(end,:)] ...
                 -  [-rho_nph(1,:); rho_nph(1:end-1,:)] )/dz/2;% Assumes rho is 0 top and bottom boundary
    
    %% Define top/bottom Neumann boundary conditions for pressure
    T = @(x) PBCT(x,time-dt/2);
    B = @(x) PBCB(x,time-dt/2);
    
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
        surf(xIn,zIn,uNew,'Edgecolor','none')
        title(['Computed u at time = ' num2str(time)])
        view([0 90]); axis tight;

        subplot(2,2,2)
        surf(xIn,zIn,wNew,'Edgecolor','none')
        title(['Computed w at time = ' num2str(time)])
        view([0 90]); axis tight;
        
        subplot(2,2,3)
        surf(xIn,zIn,s_nph,'Edgecolor','none')
        title(['Computed rho at time = ' num2str(time+dt/2)])
        view([0 90]); axis tight;
        
        subplot(2,2,4)
        surf(x,z,P,'Edgecolor','none')
        title(['Computed P at time = ' num2str(time-dt/2)])
        view([0 90]); axis tight;
        drawnow
        pause(.1)
    end
end
