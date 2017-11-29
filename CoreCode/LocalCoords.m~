function [normal,LCoords,theta] = LocalCoords(info,phi,x,y,projection)

    normal = zeros(2,size(projection,2));
    LCoords = zeros(5,size(projection,2));
    theta = zeros(1,size(projection,2));
    hx = x(1,2) - x(1,1);
    hy = y(2,1) - y(1,1);
    
    for k = 1:size(projection,2)
        %% Extract particular projection
        proj = projection(:,k);
        
        %% Calculate normal vector
        % Determine bounding nodes of projection
        bnd_x = zeros(1,2); bnd_y = zeros(1,2);
        [tmp,bnd_x(1)] = min(abs(x(1,:)-proj(1)));% identify closest grid node in x
        
        if bnd_x(1) == 1% Identify next closest
            [tmp,bnd_x(2)] = min(abs([x(1,end) x(1,bnd_x(1)+1)]-proj(1)));
            if bnd_x(2) == 1
                bnd_x(2) = length(x(1,:));
            else 
                bnd_x(2) = bnd_x(1) + 1;
            end
        elseif bnd_x(1) == length(x(1,:))
            [tmp,bnd_x(2)] = min(abs([x(1,bnd_x(1)-1) x(1,1)]-proj(1)));
            if bnd_x(2) == 1
                bnd_x(2) = bnd_x(1) - 1;
            else 
                bnd_x(2) = 1;
            end
        else
            [tmp,bnd_x(2)] = min(abs([x(1,bnd_x(1)-1) x(1,bnd_x(1)+1)]-proj(1)));
            if bnd_x(2) == 1
                bnd_x(2) = bnd_x(1) - 1;
            else 
                bnd_x(2) = bnd_x(1) + 1;
            end
        end
        
        [tmp,bnd_y(1)] = min(abs(y(:,1)-proj(2)));% identify closest grid node in y
        
        if bnd_y(1) == 1% Identify next closest
            bnd_y(2) = 2;
        elseif bnd_y(1) == length(y(:,1))
            bnd_y(2) = length(y(:,1))-1;
        else
            [tmp,bnd_y(2)] = min(abs([y(bnd_y(1)-1,1) y(bnd_y(1)+1,1)]-proj(2)));
            if bnd_y(2) == 1
                bnd_y(2) = bnd_y(1) - 1;
            else 
                bnd_y(2) = bnd_y(1) + 1;
            end
        end
        
        bnd_x = sort(bnd_x);
        bnd_y = sort(bnd_y);
        [bnd_y,bnd_x] = meshgrid(bnd_y,bnd_x);
        
        bnd_x = bnd_x(:);% column indexes for bnding box
        bnd_y = bnd_y(:);% row indexes for bnding box
        
        % Determine partial derivatives at bounding nodes
        for l = 1:length(bnd_x)
            if bnd_x(l) == 1
                dx_bnd(l) = (phi(bnd_y(l),bnd_x(l)+1)-phi(bnd_y(l),end))/(2*hx);% x-deriv
            elseif bnd_x(l) == length(x(1,:))
                dx_bnd(l) = (phi(bnd_y(l),1)-phi(bnd_y(l),bnd_x(l)-1))/(2*hx);% x-deriv
            else
                dx_bnd(l) = (phi(bnd_y(l),bnd_x(l)+1)-phi(bnd_y(l),bnd_x(l)-1))/(2*hx);% x-deriv
            end
            dy_bnd(l) = (phi(bnd_y(l)+1,bnd_x(l))-phi(bnd_y(l)-1,bnd_x(l)))/(2*hy);% y-deriv
        end
        
        %% Bilinear Interpolant
        [idx] = sub2ind(size(phi),bnd_y,bnd_x);
        dx_interp = bilinear_interp(x(idx),y(idx),dx_bnd,hx,hy);
        dy_interp = bilinear_interp(x(idx),y(idx),dy_bnd,hx,hy);
        
        %% Calculate normal
        normal(:,k) = [dx_interp(proj(1),proj(2)),dy_interp(proj(1),proj(2))];
        normal(:,k) = normal(:,k)/norm(normal(:,k));

        %% Calculate theta(angle between normal and x-axis)
        theta(k) = atan(normal(2,k)/normal(1,k));
        
        %% Calculate local coordinates for "+" and "-" nodes
        LCoords(:,k) = (info(:,6,k) - proj(1))*cos(theta(k)) + ...
            (info(:,7,k) - proj(2))*sin(theta(k));
    end
    
    %plus_Lnodes_xi = (x_Lstencil(plus_Lnodes)-proj(1))*cos(theta) + ...
    %    (y_Lstencil(plus_Lnodes)-proj(2))*sin(theta);% For laplacian stencil
    %minus_Lnodes_xi = (x_Lstencil(minus_Lnodes)-proj(1))*cos(theta) + ...
    %    (y_Lstencil(minus_Lnodes)-proj(2))*sin(theta);% For laplacian stencil
    %plus_Dnodes_xi = (x_Dstencil(plus_Dnodes)-proj(1))*cos(theta) + ...
    %    (y_Dstencil(plus_Dnodes)-proj(2))*sin(theta);% For drhodz stencil
    %minus_Dnodes_xi = (x_Dstencil(minus_Dnodes)-proj(1))*cos(theta) + ...
    %    (y_Dstencil(minus_Dnodes)-proj(2))*sin(theta);% For drhodz stencil
    
end