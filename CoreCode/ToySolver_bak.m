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
%        sv = boolian save parameter
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

function [x,z,dx,dz,time,uNew,wNew,P,s_nph,rho_nph] = ToySolver(input,Lx,Lz,m,n,tFinal,dt,vis,sv)

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

%% Define exact solution through input file
input()

%% Initial condition for velocity
uOld = uInit(xIn,zIn);
wOld = wInit(xIn,zIn);

%% Obtain s_n+1/2 using one step of Heun's method
s0 = sInit(xIn,zIn);% define initial s
rho0 = rho(xIn,zIn,s0);% define initial rho

%% Calculate rho_z as source for initial pressure
% Assumes rho is 0 top and bottom boundary
%dzrho0 = ( [rho0(2:end,:); -rho0(end,:)] ...
%                 -  [-rho0(1,:); rho0(1:end-1,:)] )/dz/2;

% Corrects so that rho is 0 top and bottom boundary adds non-zero part back
dzrho0 = ([rho0(2,:)+cm*bbeta*(zIn(2,:)-h); rho0(3:end,:); -(rho0(end,:)+cp*bbeta*(zIn(end,:)-h))] - [-(rho0(1,:)+cm*bbeta*(zIn(1,:)-h)); rho0(1:end-2,:); rho0(end-1,:)+cp*bbeta*(zIn(end-1,:)-h)])/dz/2;
dzrho0(1,:) = dzrho0(1,:) - bbeta*cm;
dzrho0(end,:) = dzrho0(end,:) - bbeta*cp;

%% Define top/bottom Neumann boundary conditions for pressure
T = @(x) PBCT(x,0);
B = @(x) PBCB(x,0);

%% Solve for Pressure at t = 0
%P = Poisson(-g*dzrho0,x,z,B,T);
P = PoissonIIM_V1(-g*dzrho0,x,z,B,T,s0);

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

    disp(['Time: ' num2str(time) '; Step: ' num2str(n) '/' num2str(N)])

    %% Compute drho/dz at n+1/2
    
    % Assumes rho is periodic
    %dzrho_nph = ([rho_nph(2:end,:); rho_nph(1,:)] - [rho_nph(end,:); ...
    %                    rho_nph(1:end-1,:)])/dz/2;
    
    % Assumes rho is 0 top and bottom boundary
    %dzrho_nph = ([rho_nph(2:end,:); -rho_nph(end,:)] ...
    %             -  [-rho_nph(1,:); rho_nph(1:end-1,:)] )/dz/2;

    % Corrects so that rho is 0 top and bottom boundary adds non-zero part back
    dzrho_nph = ([rho_nph(2,:)+cm*bbeta*(zIn(2,:)-h); ...
                  rho_nph(3:end,:); -(rho_nph(end,:)+cp*bbeta*(zIn(end,:)-h))]...
                 - [-(rho_nph(1,:)+cm*bbeta*(zIn(1,:)-h)); rho_nph(1:end-2,:);...
                    rho_nph(end-1,:)+cp*bbeta*(zIn(end-1,:)-h)])/dz/2;
    dzrho_nph(1,:) = dzrho_nph(1,:) - bbeta*cm;
    dzrho_nph(end,:) = dzrho_nph(end,:) - bbeta*cp;

    %% Solve for Pressure at n+1/2
    %P = Poisson(-g*dzrho_nph,x,z,B,T);
    P = PoissonIIM_V1(-g*dzrho_nph,x,z,B,T,s_nph);
    
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

    if n == 1
        uMax = max(uNew(:));
        uMin = min(uNew(:));
        wMax = max(wNew(:));
        wMin = min(wNew(:));
    end

    uMaxOut = max(uNew(:));
    disp(['uMax: ' num2str(uMaxOut)])
    uError = uNew - uExact(xIn,zIn,time);
    disp(['uError: ' num2str(max(uError(:)))])
    wError = wNew - wExact(xIn,zIn,time);
    disp(['wError: ' num2str(max(wError(:)))])
    SError = s_nph - SExact(xIn,zIn,time+dt/2);
    disp(['SError: ' num2str(max(SError(:)))])
    %% Correct for unknown constant
    const = PExact(x(5,5),z(5,5),time-dt/2) - P(5,5);
    P = P + const;
    PError = P(2:end-1,:) - PExact(xIn,zIn,time-dt/2);
    disp(['PError: ' num2str(max(PError(:)))])

    if vis
        %% Plot velocity and error
        fig = figure(1);
        
        subplot(2,2,1)
        pcolor(xIn,zIn,uNew)
        title(['Computed u at time = ' num2str(time)])
        shading flat; axis equal; colorbar; axis([0 Lx 0 Lz]);
        caxis([uMin,uMax])

        subplot(2,2,2)
        pcolor(xIn,zIn,wNew)
        title(['Computed w at time = ' num2str(time)])
        shading flat; axis equal; colorbar; axis([0 Lx 0 Lz]);
        caxis([uMin,uMax])

        subplot(2,2,3)
        pcolor(xIn,zIn,s_nph)
        title(['Computed s at time = ' num2str(time+dt/2)])
        shading flat; axis equal; colorbar; axis([0 Lx 0 Lz]);
        hold on;
        contour(xIn,zIn,s_nph,[0 0],'r-')
        hold off;
        
        subplot(2,2,4)
        pcolor(x,z,P)
        title(['Computed P at time = ' num2str(time-dt/2)])
        shading flat; axis equal; colorbar; axis([0 Lx 0 Lz]);

        drawnow
        
        if 0
            figure(2)
            subplot(2,2,1)
            surf(xIn,zIn,uError,'edgecolor','none')
            xlabel('x')
            ylabel('z')
            title('Error in u')
            
            subplot(2,2,2)
            surf(xIn,zIn,wError,'edgecolor','none')
            xlabel('x')
            ylabel('z')
            title('Error in w')
            
            subplot(2,2,3)
            surf(xIn,zIn,SError,'edgecolor','none')
            xlabel('x')
            ylabel('z')
            title('Error in S')
            
            subplot(2,2,4)
            surf(xIn,zIn,PError,'edgecolor','none')
            xlabel('x')
            ylabel('z')
            title('Error in P')
            drawnow;
            %pause
        end
    end
    
    if (~mod(n,round(N/20)) | n==1) & sv
        str = ['AsymSol_t' num2str(time) '.jpg'];
        saveas(fig,str);
        disp(['Saved at time t = ' num2str(time)])
    end
    
end
