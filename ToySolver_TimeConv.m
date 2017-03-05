%% Main driving code for WaveTank project

clear all; close all;

%% Define input file
input = @input_PerXNeuY;
input();

%% Define spatial mesh parameters
Lx = 2*pi;
Lz = 2*pi;
m = 2^(9);
n = 2^(9)+2;

%% for loop for temporal convergence
for i = 1:4
    
    %% Create temporal computing mesh
    tFinal = 1.0;
    dt = tFinal/2^(0+i);
    
    %% Display mesh info
    disp(['Mesh ' num2str(i) ', dt = ' num2str(dt)])

    %% Turn plotting on/off
    vis = true;
    
    %% Call solver
    [x,z,dx,dz,time,u,w,P,rho] = ToySolver(input,Lx,Lz,m,n,tFinal,dt,vis);
    
    %% Compute point-wise errors
    uError = u - uExact(x(2:end-1,:),z(2:end-1,:),time);
    wError = w - wExact(x(2:end-1,:),z(2:end-1,:),time);
    rhoError = rho - rhoExact(x(2:end-1,:),z(2:end-1,:),time+dt/2);
    PError = P - PExact(x,z,time-dt/2);    

    %% Compute \infty-norm errors
    uInfError(i) = max(abs(uError(:)));
    wInfError(i) = max(abs(wError(:)));
    rhoInfError(i) = max(abs(rhoError(:)));
    TimeStep(i) = dt;
    disp(['u Inf-norm error: ' num2str(uInfError(i))])
end

%% Plot convergence results
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