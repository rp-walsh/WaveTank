%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script tests the elliptic solve for the waterless cloud example.
%That is, we solve the problem L[P] = \rho_z with periodic BCs in x and
%Neumann in y
%with periodic B.C.'s in x and neumann in y.
%Note that this solver assumes a staggered grid in y.
%Parameters:
%           sol=solution vector
%           rhs=rhs vector for L(u)=f
%           m,n=number of points in x,y respectively
%           Lx,Ly=length of domain in x,y respectively
%           B=function call for bottom BC
%           T=function call for top BC
%
% Ray Walsh -- 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

% Add core code directory
addpath(genpath('CoreCode'));

F = @(x,z) exp(sin(2*pi*x) + z);
%source = @(x,z) 2*pi^2*exp(sin(2*pi*x) + z).*(-2*sin(2*pi*x) + cos(4*pi*x) + 1);
source = @(x,z) 4*pi^2*exp(sin(2*pi*x) + z).*(cos(2*pi*x).^2 - sin(2*pi*x) + 1/(4*pi^2));
BCT = @(x) F(x,1);
%BCT = @(x) zeros(size(x));
BCB = @(x) F(x,0);
%BCB = @(x) zeros(size(x));

for j = 1:5%for loop for convergence study

    disp(['Mesh: ' num2str(j)])
    
    %% Set up finite difference grid
    Lx = 1; Lz=1;
    if j == 1
        m = 2^(j+4)-2; n=2^(j+4);% was 4
        %m = 4; n = 6;
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

    sol = Poisson(source(x(2:end-1,:),z(2:end-1,:)),x,z,BCB,BCT);

    %% Correct for unknown constant
    const = F(x(1),z(1)) - sol(1,1);
    sol = sol + const;

    error(j) = max(abs(sol(:) - F(x(:),z(:))));
    h(j) = hx;
    %(sol(end,:) - sol(end-1,:))/hz
    %AA = BCT(x(end,:))
    %sol(end,end)

    %FD = F(x,z);
    %tmp = (FD(2:end-1,3:end) - 2*FD(2:end-1,2:end-1) + FD(2:end-1,1:end-2))/hx^2 ...
    %      + (FD(3:end,2:end-1) - 2*FD(2:end-1,2:end-1) + FD(1:end-2,2:end-1))/hz^2;

    %error = tmp - source(x(2:end-1,2:end-1),z(2:end-1,2:end-1));

    %figure(1000)
    %surf(x(2:end-1,2:end-1),z(2:end-1,2:end-1),tmp,'edgecolor','none')
    %figure(1001)
    %surf(x,z,source(x,z),'edgecolor','none')
    %figure(1002)
    %surf(x(2:end-1,2:end-1),z(2:end-1,2:end-1),error,'edgecolor','none')

end

figure
surf(x,z,sol,'edgecolor','none')
xlabel('x')
ylabel('y')
title('solution')

figure
surf(x,z,F(x,z),'edgecolor','none')
xlabel('x')
ylabel('y')
title('Exact')

figure
plot(log10(h),log10(error),'o')
P = polyfit(log10(h),log10(error),1);
hold on
plot(log10(h),P(1)*log10(h)+P(2));
legend('error',[num2str(P(1)) 'x + (' num2str(P(2)) ')'])
