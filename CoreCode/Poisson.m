function [sol]=Poisson(f,x,y,B,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function provides the solution to the poisson problem 
%L(sol)=rhs with periodic B.C.'s in x and neumann in y.
%Note that this solver assumes a staggered grid in y.
%Parameters:
%           sol=solution vector
%           f=rhs evaluated at all nodes except staggered(matrix form)
%           m,n=number of points in x,y respectively
%           Lx,Ly=length of domain in x,y respectively
%           B=function call for bottom BC
%           T=function call for top BC
%
% Ray Walsh -- 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract key info from inputs
m = size(x,2);
n = size(x,1);
hx = abs(x(1,1)-x(1,2));
hy = abs(y(1,1)-y(2,1));

%% Set up index matrix
NumUnk = m*n; %total number of unknowns
G = zeros(size(x));
G(1:NumUnk) = 1:NumUnk; %index matrix

%% Generate standard discrete laplacian
% construct d^2/dx^2 matrix
e = ones(m,1);
I = sparse(eye(n));
A = spdiags([e -2*e e],[-1 0 1],m,m);%1D d^2/dx^2 matrix
A(1,end) = 1; A(end,1) = 1;%correct for periodic BC's in x
dx = kron(A,I)/hx^2;%2D d^2/dx^2 matrix corrected for periodic BC's
%dx(G(1,:),:)=zeros(size(dx(G(1,:),:)));%Remove entries for points outside
dx(G(1,:),:)=sparse(size(dx(G(1,:),:),1),size(dx(G(1,:),:),2));%Remove entries for points outside
%dx(G(end,:),:)=zeros(size(dx(G(end,:),:)));%domain correcting for neumann BC's in y
dx(G(end,:),:)=sparse(size(dx(G(end,:),:),1),size(dx(G(end,:),:),2));%domain correcting for neumann BC's in y
% construct d^2/dy^2 matrix
e = ones(n,1);
I = sparse(eye(m));
A = spdiags([e -2*e e],[-1 0 1],n,n)/hy^2;%1D d^2/dy^2 matrix
A(1,1) = -1/hy; A(1,2) = 1/hy;%Correct for neumann BC's
A(end,end-1) = -1/hy; A(end,end) = 1/hy;%Correct for neumann BC's
dy = kron(I,A);%%2D d^2/dx^2 matrix

% add dx and dy for laplacian
Lapl = dx + dy;%discrete laplacian (dirichlet BC's)
Lapl(end,end) = 0;%fix one value of sol = 0 for uniquness

%% Set up rhs vector 
BCx = x(1,:);%Determine x dependent boundary nodes
BCB = B(BCx);%Define bottom boundary condition
BCT = T(BCx);%Define top boundary condition
BCT(end) = 0;%fix one value of sol=0 for uniquness
rhs = zeros(size(x));%initialize rhs

%evalute f excluding top and bottom nodes for neumann condition
%rhs(G(2:end-1,:)) = f(x(G(2:end-1,:)),y(G(2:end-1,:)));
rhs(G(2:end-1,:)) = f;
rhs(1,:) = BCB;%add bottom BC's
rhs(end,:) = BCT;%add top BC's
rhs=reshape(rhs,NumUnk,1);%reshape tp vector for linear solve

%%Linear solve for problem solution
sol = Lapl\rhs;

%%Reshape for visualization
sol = reshape(sol,n,m);
