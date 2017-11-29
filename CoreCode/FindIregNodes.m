%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function determines the irregular nodes (5-pt stencil crosses
% interface) where the interface is described by the 0 level set of some
% function phi(x,y).
% Input:
%       phi = gridded level set function defined on all nodes in grid
%       typ = string indicator for which stencil we're interested in
% Output:
%        IregNodes = linear indexes of irregular nodes
% Ray Walsh -- 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IregNodes] = FindIregNodes(phi,typ)

    idtfy = phi >= 0;%Identify points in Omega+ 

    switch lower(typ)
      case 'lapl'
        
        shiftL = circshift(idtfy,[0,-1]);%Shift points included in 5-point stencil
        shiftR = circshift(idtfy,[0,1]);%Shift points included in 5-point stencil
        shiftU = circshift(idtfy,[-1,0]);%Shift points included in 5-point stencil
        shiftD = circshift(idtfy,[1,0]);%Shift points included in 5-point stencil
        indctr = [idtfy(1,:) + idtfy(2,:) + 5; idtfy(2:end-1,:) + shiftL(2:end-1,:)...
                  + shiftR(2:end-1,:) + shiftU(2:end-1,:) + shiftD(2:end-1,:);...
                  idtfy(end,:) + idtfy(end-1,:) + 5];%Create matrix summing values in stencil
        IregNodes = find(indctr~=0 & indctr~=5 & indctr~=7);%Determine irregular nodes
        clear idtfy indctr shiftL shiftR shiftU shiftD;
        
      case 'dx'
        
        shiftL = circshift(idtfy,[0,-1]);%Shift points included in 5-point stencil
        shiftR = circshift(idtfy,[0,1]);%Shift points included in 5-point stencil
        indctr = idtfy + shiftL + shiftR;%Create matrix summing values in stencil
        IregNodes = find(indctr~=0 & indctr~=3);%Determine irregular nodes
        clear idtfy indctr shiftL shiftR;
    
      case 'dz'
        
        shiftU = circshift(idtfy,[-1,0]);%Shift points included in 5-point stencil
        shiftD = circshift(idtfy,[1,0]);%Shift points included in 5-point stencil
        indctr = idtfy + shiftU + shiftD;%Create matrix summing values in stencil
        indctr(1,:) = zeros(size(indctr(1,:)));%bottom nodes can't be irregular
        indctr(end,:) = zeros(size(indctr(end,:)));%top nodes can't be irregular
        IregNodes = find(indctr~=0 & indctr~=3);%Determine irregular nodes
        clear idtfy indctr shiftU shiftD;
        
    end
end