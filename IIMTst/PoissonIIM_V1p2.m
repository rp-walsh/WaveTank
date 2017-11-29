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
function [sol,varargout]=PoissonIIM(f,x,y,B,T,phi,gamma,varargin)
IIM = true;

if ~isempty(varargin)
    jump = varargin{1};
end

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

%% Fix singular matrix by removing row and column of 0ed entry
%Lapl(end,:) = [];
%Lapl(:,end) = [];

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
    %proj = OrthogProjInfo(info,phi,x,y);
    proj = x(IregNodes)';
    proj(2,:) = gamma(proj(1,:));
    %figure(555)
    %mesh(x,y,zeros(size(x)));
    %hold on
    %plot(x(IregNodes),y(IregNodes),'ko');
    %plot(x(1,:),gamma(x(1,:)),'r')
    %hold off
    %view([0,90]);
    %pause
    
    %% Compute normal vector at each projection
    [normal,LCoords,theta] = LocalCoords(info,phi,x,y,proj);
    
    %% Extrapolate for jump information
    IregNodesInterp = FindIregNodes(phi(2:end-1,:),'lapl');
    MinusNodes = setdiff(find(phi(2:end-1,:)<=0),IregNodesInterp);
    PlusNodes = setdiff(find(phi(2:end-1,:)>0),IregNodesInterp);
    xInterp = x(2:end-1,:);
    yInterp = y(2:end-1,:);
    F0 = scatteredInterpolant(xInterp(MinusNodes),yInterp(MinusNodes),f(MinusNodes),'natural','nearest');
    F1 = scatteredInterpolant(xInterp(PlusNodes),yInterp(PlusNodes),f(PlusNodes),'natural','nearest');

    if isempty(varargin)
        varargout{1} = @(x,z) F1(x,z) - F0(x,z);
    end
        
    %% Calulate corrections
    coeff = [-2/hx^2-2/hy^2; 1/hy^2; 1/hx^2; 1/hy^2; 1/hx^2];
    Dcoeff = [0; 1/(2*hy); -1/(2*hy)];
    CorVec = zeros(size(rhs));
    for k = 1:size(info,3)
        plus_Lnodes_xi = LCoords(logical(info(:,4,k)),k).';
        a8 = (plus_Lnodes_xi.^2/2)*coeff(logical(info(:,4,k)));
        DNodesInfo = info([true true false true false],:,k);
        LCoordsDNodes = LCoords([true true false true false],k);
        plus_Dnodes_xi = LCoordsDNodes(logical(DNodesInfo(:,4))).';
        Dcorr = plus_Dnodes_xi*Dcoeff(logical(DNodesInfo(:,4)));
        side = info(1,5,k);
        if isempty(varargin)
            if side == 0
                correction = a8*(F1(proj(1,k),proj(2,k))-F0(proj(1,k),proj(2,k))) -...
                    (Dcorr/sin(theta(k)))*(F1(proj(1,k),proj(2,k))- ...
                                           F0(proj(1,k),proj(2,k)));
            else
                correction = a8*(F0(proj(1,k),proj(2,k))-F1(proj(1,k),proj(2,k))) -...
                    (Dcorr/sin(theta(k)))*(F0(proj(1,k),proj(2,k))- ...
                                           F1(proj(1,k),proj(2,k)));
            end
        else
            if side == 0
                correction = a8*(jump(proj(1,k),proj(2,k))) -...
                    (Dcorr/sin(theta(k)))*(jump(proj(1,k),proj(2,k)));
            else
                correction = a8*(-jump(proj(1,k),proj(2,k))) -...
                    (Dcorr/sin(theta(k)))*(-jump(proj(1,k),proj(2,k)));
            end
        end
        %rhs(info(1,1,k)) = rhs(info(1,1,k)) + correction;%For problem 1
        CorVec(info(1,1,k)) = correction;%For problem 1
    end
    %solveabilityRHS = hx*hy*sum(sum(rhs(2:end-1,:))) - hx*sum(rhs(end,:)) + hx*sum(rhs(1,:))
    %solveabilityCorVec = hx*hy*sum(sum(CorVec(2:end-1,:))) - hx*sum(CorVec(end,:)) + hx*sum(CorVec(1,:))
    CorVecMean = mean(nonzeros(CorVec(:)));
    CorVec(abs(CorVec)>0) = CorVec(abs(CorVec)>0) - CorVecMean;
    rhs = rhs + CorVec;%For problem 1
    CorVec=reshape(CorVec,NumUnk,1);%reshape tp vector for linear solve
    %soltmp = Lapl\CorVec;
    %soltmp = reshape(soltmp,n,m);
    
    %figure(101)
    %plot(y(:,1),soltmp(:,round(size(x,2)/2)),'b.-')
    %title(['z-slice of Numerical Solution for Corrections'])
    %axis tight
end

%% Reshape rhs for linear solve
%solveability = hx*hy*sum(sum(rhs(2:end-1,:))) - hx*sum(rhs(end,:)) + hx*sum(rhs(1,:))
%rhs(end) = rhs(end) + solveability/hx;
%solveability = hx*hy*sum(sum(rhs(2:end-1,:))) - hx*sum(rhs(end,:)) + hx*sum(rhs(1,:))
rhs=reshape(rhs,NumUnk,1);%reshape tp vector for linear solve

%%Linear solve for problem solution
sol = Lapl\rhs;
%sol = [sol; 0];

%%Reshape for visualization
sol = reshape(sol,n,m);