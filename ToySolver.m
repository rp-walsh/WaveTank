%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This solver uses the split-step scheme 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,z,dx,dz,uNew,wNew,P,rho_nph] = ToySolver(input,Lx,Lz,m,n,tFinal,dt,vis)

%% Define exact solution through input file
%input_PerXNeuY2;
input()

%% Define initial conditions
u0 = @(x,z) uExact(x,z,0);
w0 = @(x,z) wExact(x,z,0);

%% Create spatial computing mesh
dx = Lx/m;% corrected for periodic grid
dz = Lz/(n-2);% corrected for staggered grid
x = 0:dx:Lx-dx;
z = 0-dz/2:dz:Lz+dz/2;
[x,z] = meshgrid(x,z);
xIn = x(2:end-1,:);
zIn = z(2:end-1,:);

%% loop for temporal convergence
%for i = 1:4

    %% Initial condition for velocity
    uOld = u0(xIn,zIn);
    wOld = w0(xIn,zIn);
    
    %% Create temporal computing mesh
    N = round(tFinal/dt);
    
    %% Obtain rho_n+1/2 (use rho(0) when no exact solution)
    rho_nph = rhoExact(xIn,zIn,dt/2);

    %% Begin time-stepping
    for n=1:N
        
        time = n*dt;
        
        %% Compute drho/dz at n+1/2
        %dzrho_nph = ([rho_nph(2:end,:); rho_nph(1,:)] - [rho_nph(end,:); ...
        %                    rho_nph(1:end-1,:)])/dz/2;% Assumes rho is periodic
        dzrho_nph = ([rho_nph(2:end,:); -rho_nph(end,:)] ...
                  -  [-rho_nph(1,:); rho_nph(1:end-1,:)] )/dz/2;% Assumes rho is 0 top and bottom boundary

        %% Define top/bottom Neumann boundary conditions for pressure
        T = @(x) PExact_z(x,2*pi,time-dt/2);
        B = @(x) PExact_z(x,0,time-dt/2);
        
        %% Solve for Pressure at n+1/2
        P = Poisson(-g*dzrho_nph,x,z,B,T);
        
        %% Correct for unknown constant
        PExactDisc = PExact(x,z,0);
        const = P(10,10) - PExactDisc(10,10);
        P = P - const;
        
        %% Compute grad P
        P_z = (P(3:end,:) - P(1:end-2,:))/2/dz;
        P_x = ([P(2:end-1,2:end) P(2:end-1,1)] - [P(2:end-1,end) P(2:end-1,1:end-1)])/2/dx;
        
        %% Step velocity forward in time
        uNew = uOld - dt*(P_x);
        wNew = wOld - dt*(P_z + g*rho_nph);
        uOld = uNew;
        wOld = wNew;
        
        %% Step density forward in time
        rho_nph = rho_nph + dt*wNew;
        
        %% Compute point-wise errors
        uError = uNew - uExact(xIn,zIn,time);
        wError = wNew - wExact(xIn,zIn,time);
        rhoError = rho_nph - rhoExact(xIn,zIn,time+dt/2);

        if vis
            %% Plot velocity and error
            figure(1)
            surf(xIn,zIn,rho_nph,'Edgecolor','none')
            title(['Computed rho at time = ' num2str(time+1/2)])
            view([0 90]);
            figure(2)
            surf(x,z,P,'Edgecolor','none')
            title(['Computed P at time = ' num2str(time+1/2)])
            view([0 90]);
            drawnow
            pause(.1)
        end
    end
    
    %% Compute Inf-norm errors
    %uInfError(i) = max(abs(uError(:)));
    %wInfError(i) = max(abs(wError(:)));
    %rhoInfError(i) = max(abs(rhoError(:)));
    %TimeStep(i) = dt;

    %disp(['Inf-norm error: ' num2str(uInfError(i))])

    %end

if 0
    figure
    plot(log(TimeStep),log(uInfError),'o')
    pfit = polyfit(log(TimeStep),log(uInfError),1);
    hold on
    plot(log(TimeStep),pfit(1)*log(TimeStep) + pfit(2));
    legend('Numerical Convergence',[num2str(pfit(1)) ' + (' num2str(pfit(2)) ')'])
    xlabel('$\log(\Delta t)$','interpreter','latex')
    ylabel('$\log(||u-u_h||_\infty)$','interpreter','latex')
    title('Temporal Convergence of u')
    
    figure
    plot(log(TimeStep),log(wInfError),'o')
    pfit = polyfit(log(TimeStep),log(wInfError),1);
    hold on
    plot(log(TimeStep),pfit(1)*log(TimeStep) + pfit(2));
    legend('Numerical Convergence',[num2str(pfit(1)) ' + (' num2str(pfit(2)) ')'])
    xlabel('$\log(\Delta t)$','interpreter','latex')
    ylabel('$\log(||w-w_h||_\infty)$','interpreter','latex')
    title('Temporal Convergence of w')
    
    figure
    plot(log(TimeStep),log(rhoInfError),'o')
    pfit = polyfit(log(TimeStep),log(rhoInfError),1);
    hold on
    plot(log(TimeStep),pfit(1)*log(TimeStep) + pfit(2));
    legend('Numerical Convergence',[num2str(pfit(1)) ' + (' num2str(pfit(2)) ')'])
    xlabel('$\log(\Delta t)$','interpreter','latex')
    ylabel('$\log(||\rho-\rho_h||_\infty)$','interpreter','latex')
    title('Temporal Convergence of \rho')
end