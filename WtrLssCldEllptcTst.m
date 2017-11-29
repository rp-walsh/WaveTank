%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tests the elliptic solve for the waterless cloud example.
% That is, we solve the problem L[P] = \rho_z with periodic BCs in x and
% Neumann in y
% with periodic B.C.'s in x and neumann in y.
% Note that this solver assumes a staggered grid in y.
% Parameters:
%            sol=solution vector
%            rhs=rhs vector for L(u)=f
%            m,n=number of points in x,y respectively
%            Lx,Ly=length of domain in x,y respectively
%            B=function call for bottom BC
%            T=function call for top BC
%
% Ray Walsh -- 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

%% Add subdirectories
addpath(genpath('Inputs'));
addpath(genpath('CoreCode'));

%% Define input file
input = @input_WtrLssCld;
input();

for j = 1:4%for loop for convergence study (was 5)

    disp(['Mesh: ' num2str(j)])
    
    %% Set up finite difference grid
    Lx = 1; Lz=1;
    if j == 1
        m = 2^(j+2)-2; n=2^(j+2);% was 4(2)
        hx = Lx/m; %define stepsize in x (periodic domain)
        hz = Lz/(n-2); %define stepsize in y (staggered)
        [x z] = meshgrid(0:hx:Lx-hx,-hz/2:hz:Lz+hz/2);
    else
        hx = hx/2; %define stepsize in x (periodic domain)
        hz = hz/2; %define stepsize in y (staggered)
        [x z] = meshgrid(0:hx:Lx-hx,-hz/2:hz:Lz+hz/2);
        m = size(x,2);
        n = size(x,1);
    end

    %% Differentiate RHS
    rho_disc = rhoExact(x,z,0);
    drho_discdz = (rho_disc(3:end,:) - rho_disc(1:end-2,:))/(2*hz);

    %% Test differentiation 
    rho0 = rho_disc(2:end-1,:);
    xIn= x(2:end-1,:);
    zIn= z(2:end-1,:);
    dzrho0 = ( [rho0(2,:)+cm*bbeta*(zIn(2,:)-h); rho0(3:end,:); -(rho0(end,:)+cp*bbeta*(zIn(end,:)-h))]...
    -  [-(rho0(1,:)+cm*bbeta*(zIn(1,:)-h)); rho0(1:end-2,:); rho0(end-1,:)+cp*bbeta*(zIn(end-1,:)-h)] )/hz/2;
    dzrho0(1,:) = dzrho0(1,:) - bbeta*cm;
    dzrho0(end,:) = dzrho0(end,:) - bbeta*cp;

    %% Calculate laplacian of pressure
    PDisc = PExact(x,z,0);
    LaplP = (PDisc(3:end,2:end-1) + PDisc(2:end-1,3:end) - 4*PDisc(2:end-1,2:end-1) ...
             + PDisc(1:end-2,2:end-1) + PDisc(2:end-1,1:end-2))/hx^2;
     
    %sol = Poisson(-g*drho_discdz,x,z,@(x) PBCB(x,0),@(x) PBCT(x,0));
    sol = Poisson(-g*dzrho0,x,z,@(x) PBCB(x,0),@(x) PBCT(x,0));

    %% Correct for unknown constant
    const = PDisc(5,5) - sol(5,5);
    sol = sol + const;

    %% Calculate error
    PoissonErr(j) = max(abs(sol(:) - PDisc(:)));
    dx(j) = hx;

    %LaplSol = (sol(3:end,2:end-1) + sol(2:end-1,3:end) - 4*sol(2:end-1,2:end-1) ...
    %         + sol(1:end-2,2:end-1) + sol(2:end-1,1:end-2))/hx^2;

    sol_z = (sol(3:end,:) - sol(1:end-2,:))/hz/2;
    sol_x = (sol(:,3:end) - sol(:,1:end-2))/hx/2;
    
    PzDisc = PzExact(x(2:end-1,:),z(2:end-1,:),0);
    ZDerivErr(j) = abs(max(sol_z(:) - PzDisc(:)));

    PxDisc = PxExact(x(:,2:end-1),z(:,2:end-1),0);
    XDerivErr(j) = abs(max(sol_x(:) - PxDisc(:)));
end

if 0
    figure
    surf(x(2:end-1,:),z(2:end-1,:),-g*drho_discdz,'edgecolor','none')
    xlabel('x')
    ylabel('y')
    title('\rho_z')
    
    figure
    surf(x(2:end-1,2:end-1),z(2:end-1,2:end-1),LaplP,'edgecolor','none')
    xlabel('x')
    ylabel('y')
    title('L[P]')
end

figure
surf(x,z,sol,'edgecolor','none')
xlabel('x')
ylabel('y')
title('Numerical solution')

figure
surf(x,z,PDisc,'edgecolor','none')
xlabel('x')
ylabel('y')
title('Exact solution')

figure
plot(log10(dx),log10(PoissonErr),'o')
Pfit = polyfit(log10(dx),log10(PoissonErr),1);
hold on
plot(log10(dx),Pfit(1)*log10(dx)+Pfit(2));
legend('error',[num2str(Pfit(1)) 'x + (' num2str(Pfit(2)) ')'])
title('Poisson solve convergence')

figure
plot(log10(dx),log10(ZDerivErr),'o')
Pfit = polyfit(log10(dx),log10(ZDerivErr),1);
hold on
plot(log10(dx),Pfit(1)*log10(dx)+Pfit(2));
legend('error',[num2str(Pfit(1)) 'x + (' num2str(Pfit(2)) ')'])
title('Z derivative convergence')

figure
plot(log10(dx),log10(XDerivErr),'o')
Pfit = polyfit(log10(dx),log10(XDerivErr),1);
hold on
plot(log10(dx),Pfit(1)*log10(dx)+Pfit(2));
legend('error',[num2str(Pfit(1)) 'x + (' num2str(Pfit(2)) ')'])
title('X derivative convergence')
