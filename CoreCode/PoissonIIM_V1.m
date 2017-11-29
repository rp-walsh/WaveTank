%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function provides the solution to the poisson problem 
%L(sol)=rhs with a jump discontinunity in rhs. The BCs are assumed
%periodic in x and neumann in y with [sol] = [sol_\hat{n}] = 0. This
%solver uses the propper IIM discretization. 
%Note that this solver assumes a staggered grid in y.
%Parameters:
%           sol = solution vector
%           f   = rhs evaluated at all nodes except staggered top and
%                 bottom nodes(matrix form)
%           x,y =  x,y mesh coordinate matricies
%           B   = function call for bottom neumann BC
%           T   = function call for top neumann BC
%           phi = Gridded level set function defining the interface
%
%Dependcies:
%           FindIregNodes.m
%           StencilInfo.m
%           OrthogProjInfo.m
%           LocalCoords.m
%           bilinear_interp.m
%
% Ray Walsh -- 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sol]=PoissonIIM(f,x,y,B,T,phi)
%IIM = true;
IIM = true;

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
Lapl(end,end-1) = 0;%fix one value of sol = 0 for uniquness
%Lapl(1,2) = 0;
%tmp = length(x(:));
%tmp = tmp/2;
%Lapl(tmp,:) = zeros(size(Lapl(tmp,:)));
%Lapl(tmp,tmp) = 1;

%% Set up rhs vector 
BCx = x(1,:);%Determine x dependent boundary nodes
BCB = B(BCx);%Define bottom boundary condition
BCT = T(BCx);%Define top boundary condition
BCT(end) = 0;%fix one value of sol=0 for uniquness
rhs = zeros(size(x));%initialize rhs

%% Evalute f excluding top and bottom nodes for neumann condition
%rhs(G(2:end-1,:)) = f(x(G(2:end-1,:)),y(G(2:end-1,:)));
rhs(G(2:end-1,:)) = f;
rhs(1,:) = BCB;%add bottom BC's
rhs(end,:) = BCT;%add top BC's

if IIM
    %%CORRECT LATER!!! QUICK FIX DUE TO SIZE OVERSIGHT ON PHI
    %%WILL NEED TO ADJUST S.T. PHI AND F CAN BE SAME SIZE
    phi = [phi(1,:); phi; phi(end,:)];
    
    %%%%%%% Compute IIM correction for rhs
    IregNodes = FindIregNodes(phi,'lapl');
    
    %% Obtain required stencil information for irregular nodes
    info = StencilInfo(IregNodes,phi,x,y);
    
    %% Compute orthogonal projection for each node
    proj = OrthogProjInfo(info,phi,x,y);
    
    %% Compute normal vector at each projection
    [normal,LCoords,theta] = LocalCoords(info,phi,x,y,proj);
    
    %% Extrapolate for jump information
    MinusNodes = find(phi(2:end-1,:)<=0);
    PlusNodes = find(phi(2:end-1,:)>0);
    xInterp = x(2:end-1,:);
    yInterp = y(2:end-1,:);
    F0 = scatteredInterpolant(xInterp(MinusNodes),yInterp(MinusNodes),f(MinusNodes),'natural','nearest');
    F1 = scatteredInterpolant(xInterp(PlusNodes),yInterp(PlusNodes),f(PlusNodes),'natural','nearest');
    
    %% Calulate corrections
    coeff = [-2/hx^2-2/hy^2; 1/hy^2; 1/hx^2; 1/hy^2; 1/hx^2];
    Dcoeff = [0; 1/(2*hy); -1/(2*hy)];
    for k = 1:size(info,3)
        plus_Lnodes_xi = LCoords(logical(info(:,4,k)),k).';
        a8 = (plus_Lnodes_xi.^2/2)*coeff(logical(info(:,4,k)));
        DNodesInfo = info([true true false true false],:,k);
        LCoordsDNodes = LCoords([true true false true false],k);
        plus_Dnodes_xi = LCoordsDNodes(logical(DNodesInfo(:,4))).';
        Dcorr = plus_Dnodes_xi*Dcoeff(logical(DNodesInfo(:,4)));
        side = info(1,5,k);
        if side == 0
            correction = a8*(F1(proj(1,k),proj(2,k))-F0(proj(1,k),proj(2,k))) -...
                (Dcorr/sin(theta(k)))*(F1(proj(1,k),proj(2,k))- ...
                                       F0(proj(1,k),proj(2,k)));
        else
            correction = a8*(F0(proj(1,k),proj(2,k))-F1(proj(1,k),proj(2,k))) -...
                (Dcorr/sin(theta(k)))*(F0(proj(1,k),proj(2,k))- ...
                                       F1(proj(1,k),proj(2,k)));
        end
        rhs(info(1,1,k)) = rhs(info(1,1,k)) + correction;%For problem 1
        %IregNodes(k)
        %if correction~=0
        %    correction
        %end
        %pause
    end

    %figure(1)
    %surf(x(2:end-1,:),y(2:end-1,:),f,'edgecolor','none')
    %title('plot of RHS')
    %figure(2)
    %plot(y(2:end-1,round(size(x,2)/4)),f(:,round(size(x,2)/4)),'.-')
    %xlabel('z')
    %ylabel('f')
    %title('1D plot of RHS')
    %disp('pause')
    %pause
end

%% Reshape rhs for linear solve
rhs=reshape(rhs,NumUnk,1);%reshape tp vector for linear solve

%%Linear solve for problem solution
sol = Lapl\rhs;

%%Reshape for visualization
sol = reshape(sol,n,m);
