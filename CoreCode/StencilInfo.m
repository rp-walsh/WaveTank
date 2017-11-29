%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to determine which nodes belong to the 5-point
% stencil of an irregular node as well as determine which nodes lie on
% which side of the interface.
%
% Inputs:
%         IregNodes = List of N linear indexs for irregular nodes
%         phi = Gridded level set function defining the interface
%         Note:- the list in IregNodes must be compatable with phi
%                i.e. the list of linear indexes must assume the same
%                dimension as phi. 
%
% Output: info = 5x5xN arrary of the form [ M_L, M_R, M_C, op, side; ]
%                                         [ N_L, N_R, N_C, op, side; ]
%                                         [ E_L, E_R, E_C, op, side; ]
%                                         [ S_L, S_R, S_C, op, side; ]
%                                         [ W_L, W_R, W_C, op, side; ]
%
%         Note:- op is a logical variable which is true for node on the
%                opposite side of the master node.
%              - side is a logical variable which is true for nodes on
%                the positive side of the interface and false for nodes
%                on the negative side of the interface
%
% Ray Walsh -- 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [info] = StencilInfo(IregNodes,phi,x,y)

    info = zeros(5,7,length(IregNodes));
    hx = x(1,2) - x(1,1);
    hy = y(2,1) - y(1,1);
    
    for k = 1:length(IregNodes)

        %% Store master node in info
        info(1,1,k) = IregNodes(k);% Store linear index of master node
        [info(1,2,k),info(1,3,k)] = ind2sub(size(phi),info(1,1,k));% find
                                                                   % &
                                                                   % store
                                                                   % subscript
                                                                   % indicies 
                                                                   % for master node
        %% Determine side of master node
        if phi(info(1,1,k)) < 0% determine which side master node is on
            side = 0;
        else
            side = 1;
        end

        %% Find North, South nodes in discrete laplacian stencil
        info(2,2:3,k) = [info(1,2,k)+1 info(1,3,k)];% define north node
        info(2,1,k) = sub2ind(size(phi),info(2,2,k),info(2,3,k));% store linear index
        info(4,2:3,k) = [info(1,2,k)-1 info(1,3,k)];% define south node
        info(4,1,k) = sub2ind(size(phi),info(4,2,k),info(4,3,k));% store linear index
        
        %% Find East, West nodes in discrete laplacian stencil
        if info(1,3,k) == 1% Find correct 5-pt stencil for periodic BC's
            info(3,2:3,k) = [info(1,2,k), info(1,3,k)+1];% define node to east
            info(3,1,k) = sub2ind(size(phi),info(3,2,k),info(3,3,k));% store linear index
            info(5,2:3,k) = [info(1,2,k), size(phi,2)];% define node to west
            info(5,1,k) = sub2ind(size(phi),info(5,2,k),info(5,3,k));% store linear index
            % define 'Proper' x-coords for stencil. Used for local coordinates.
            info(:,6,k) = [x(info(1,1,k));x(info(2,1,k));x(info(3,1,k));x(info(4,1,k));x(info(1,1,k))-hx];
        elseif  info(1,3,k) == size(phi,2)% Find correct 5-pt stencil for periodic BC's
            info(3,2:3,k) = [info(1,2,k) 1];% define node to east
            info(3,1,k) = sub2ind(size(phi),info(3,2,k),info(3,3,k));% store linear index
            info(5,2:3,k) = [info(1,2,k) info(1,3,k)-1];% define node to west
            info(5,1,k) = sub2ind(size(phi),info(5,2,k),info(5,3,k));% store linear index
            info(:,6,k) = [x(info(1,1,k));x(info(2,1,k));x(info(1,1,k))+hx;x(info(4,1,k));x(info(5,1,k))];
            %x_Lstencil = [x(M(3)) x(L(3)) x(M(3))+hx x(M(3)) x(M(3))];% see above
        else
            info(3,2:3,k) = [info(1,2,k) info(1,3,k)+1];% define node to east
            info(3,1,k) = sub2ind(size(phi),info(3,2,k),info(3,3,k));% store linear index
            info(5,2:3,k) = [info(1,2,k) info(1,3,k)-1];% define node to west
            info(5,1,k) = sub2ind(size(phi),info(5,2,k),info(5,3,k));% store linear index
            info(:,6,k) = [x(info(1,1,k)); x(info(2,1,k)); x(info(3,1,k)); x(info(4,1,k)); x(info(5,1,k))];
            %x_Lstencil = [x(M(3)) x(L(3)) x(R(3)) x(M(3)) x(M(3))];% see above
        end

        info(:,7,k) = [y(info(1,1,k)); y(info(2,1,k)); y(info(3,1,k)); y(info(4,1,k)); y(info(5,1,k))];
        %y_Lstencil = [y(M(3)) y(M(3)) y(M(3)) y(Bot(3)) y(Top(3))];% see above
        
        %% Create stencils for drho/dz
        %x_Dstencil = [x_Lstencil(1) x_Lstencil(end-1) x_Lstencil(end)];
        %y_Dstencil = [y_Lstencil(1) y_Lstencil(end-1) y_Lstencil(end)];
        
        %% Determine nodes on opposite side of interface from master
        %Lstencil = [M(3), L(3), R(3), Bot(3), Top(3)];% Define stencil vector laplacian
        %Dstencil = [Lstencil(1) Lstencil(end-1) Lstencil(end)];% Define stencil
                                                               % vector drho/dz
        %minus_Lnodes = find(rl(x(Lstencil),y(Lstencil))*...
        %                    rl(x(Lstencil(1)),y(Lstencil(1))) >= 0);% find minus
                                                                    % nodes laplacian
        %minus_Dnodes = find(rl(x(Dstencil),y(Dstencil))*...
        %                    rl(x(Dstencil(1)),y(Dstencil(1))) >= 0,1,'last');% find
                                                                             % minus
                                                                             % nodes drho/dz
        info(:,4,k) = phi(info(1,1,k))*phi(info(:,1,k)) < 0;% find op nodes
        info(:,5,k) = phi(info(:,1,k)) > 0;
        %info(:,4,k) = plus_nodes;
        %plus_Dnodes = find(rl(x(Dstencil),y(Dstencil))*...
        %                   rl(x(Dstencil(1)),y(Dstencil(1))) < 0,1,'last');% find plus nodes
    end
    
end