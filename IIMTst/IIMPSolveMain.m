%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tests the elliptic IIM solve for eventual implimentation in
% WaveTank.
% That is, we solve the problem L[P] = \rho_z with periodic BCs in x and
% Neumann in y using the proper IIM.
% Note that this solver assumes a staggered grid in y.
% Parameters:
%            sol=solution vector
%            rhs=rhs vector for L(u)=f
%            m,n=number of points in x,y respectively
%            Lx,Ly=length of domain in x,y respectively
%            B=function call for bottom BC
%            T=function call for top BC
%
% Notes:- Differentiating rho for the RHS has been modified for the
%         masters tst case.
%
% Ray Walsh -- 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%% Add subdirectories
addpath(genpath('../Inputs'));
addpath(genpath('../CoreCode'));

%% Define input file
%input = @input_WtrLssCld;
input = @input_EllipticTstMstrs;
input();

for j = 1:4% for loop for convergence study (was 5)

    disp(['Mesh: ' num2str(j)])
    
    %% Set up finite difference grid
    Lx = 1; Lz=1;
    if j == 1
        m = 2^(j+5)-2; n=2^(j+5);% was 4
        hx = Lx/m; %define stepsize in x (periodic domain)
        hz = Lz/(n-2); %define stepsize in y (staggered)
        [x z] = meshgrid(0:hx:Lx-hx,-hz/2:hz:Lz+hz/2);
        xIn = x(2:end-1,:);
        zIn = z(2:end-1,:);
    else
        hx = hx/2; %define stepsize in x (periodic domain)
        hz = hz/2; %define stepsize in y (staggered)
        [x z] = meshgrid(0:hx:Lx-hx,-hz/2:hz:Lz+hz/2);
        m = size(x,2);
        n = size(x,1);
        xIn = x(2:end-1,:);
        zIn = z(2:end-1,:);
    end
    
    %% Differentiate RHS
    rho_disc = rhoExact(x,z,0);
    drho_discdz = (rho_disc(3:end,:) - rho_disc(1:end-2,:))/(2*hz);

    %PApprx = PoissonIIM_V1p2(-g*drho_discdz,x,z,@(x) PBCB(x,0),@(x) PBCT(x,0),rl(xIn,zIn),jump);
    [PApprx,jump] = PoissonIIM_V1p2(-g*drho_discdz,x,z,@(x) PBCB(x,0),@(x) PBCT(x,0),rl(xIn,zIn),d);
    
    %% Correct for unknown constant
    PDisc = PExact(x,z,0);
    const = PDisc(5,5) - PApprx(5,5);
    PApprx = PApprx + const;
    
    %% Calculate error
    PoissonErr = PApprx - PDisc;
    [PoissonErrMax(j)] = max(abs(PoissonErr(:)));
    dx(j) = hx;
    dz(j) = hz;

    %% Compute Gradient
    PApprx_z = (PApprx(3:end,:) - PApprx(1:end-2,:))/hz/2;
    PApprx_x = (PApprx(:,3:end) - PApprx(:,1:end-2))/hx/2;
    [PApprx_x,PApprx_z] = Gradient(PApprx,hz,rl(xIn,zIn),jump,x,z,d);

    %% Compute error on all nodes
    PzDisc = PzExact(x(2:end-1,:),z(2:end-1,:),0);
    ZDerivErr = PApprx_z - PzDisc;
    [ZDerivErrMax(j),loc] = max(abs(ZDerivErr(:)));
    PxDisc = PxExact(x(:,2:end-1),z(:,2:end-1),0);
    XDerivErr = PApprx_x - PxDisc;
    XDerivErrMax(j) = abs(max(XDerivErr(:)));

    disp(['max error in Pz: (' num2str(xIn(loc)) ',' num2str(zIn(loc)) ')'])
    
    %% Find irregular nodes in x and z
    phi = SExact(x(2:end-1,:),z(2:end-1,:),0);
    IregNodesZ = FindIregNodes(phi,'dz');
    phi = SExact(x(:,2:end-1),z(:,2:end-1),0);
    IregNodesX = FindIregNodes(phi,'dx');

    ZDerivErrMean(j) = mean(abs(ZDerivErr(IregNodesZ)));
    
    %% Compute error on regular nodes
    PApprx_z(IregNodesZ) = nan;
    ZDerivErrReg(j) = abs(max(PApprx_z(:) - PzDisc(:)));
    PApprx_x(IregNodesX) = nan;
    XDerivErrReg(j) = abs(max(PApprx_x(:) - PxDisc(:)));

    figure(101)
    surf(x,z,PoissonErr,'edgecolor','none')
    xlabel('x')
    ylabel('y')
    title('Error in Poisson solve for P')

    figure(102)
    surf(xIn,zIn,ZDerivErr,'edgecolor','none')
    xlabel('x')
    ylabel('y')
    title('Error in P_z')
    view([0 90])
    
    % drawnow
    % pause
end

figure(1)
surf(x,z,PApprx,'edgecolor','none')
xlabel('x')
ylabel('y')
title('Approximated P')

figure(2)
surf(x,z,PDisc,'edgecolor','none')
xlabel('x')
ylabel('y')
title('Exact P')

figure(3)
surf(x,z,PoissonErr,'edgecolor','none')
xlabel('x')
ylabel('y')
title('Error in Poisson solve for P')

figure(4)
plot(log10(dx),log10(PoissonErrMax),'o')
Pfit = polyfit(log10(dx),log10(PoissonErrMax),1);
hold on
plot(log10(dx),Pfit(1)*log10(dx)+Pfit(2));
hold off
legend('error',[num2str(Pfit(1)) 'x + (' num2str(Pfit(2)) ')'])
title('Convergence in Poisson solve for P')
xlabel('log_{10}(h)')
ylabel('log_{10}(||P-P_h||_\infty)')

figure(5)
plot(log10(dz),log10(ZDerivErrMax),'o')
Pfit = polyfit(log10(dz),log10(ZDerivErrMax),1);
hold on
plot(log10(dz),Pfit(1)*log10(dz)+Pfit(2));
hold off
legend('error',[num2str(Pfit(1)) 'x + (' num2str(Pfit(2)) ')'])
title('Convergence in P_z')
xlabel('log_{10}(h)')
ylabel('log_{10}(||Pz-Pz_h||_\infty)')

figure(6)
plot(log10(dz),log10(ZDerivErrMean),'o')
Pfit = polyfit(log10(dz),log10(ZDerivErrMean),1);
hold on
plot(log10(dz),Pfit(1)*log10(dz)+Pfit(2));
hold off
legend('error',[num2str(Pfit(1)) 'x + (' num2str(Pfit(2)) ')'])
title('Convergence in P_z')
xlabel('log_{10}(h)')
ylabel('log_{10}(||Pz-Pz_h||_{mean})')
