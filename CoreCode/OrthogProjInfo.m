function proj = OrthogProjInfo(info,phi,x,z)
    
    %% Extract mesh spacing 
    hx = x(1,2) - x(1,1);
    hz = z(2,1) - z(1,1);

    %% Preallocate proj output
    proj = zeros(2,size(info,3));
    
    %% Loop through irregular nodes
    for k = 1:length(proj)
        
        %% Extract required information
        node = info(1,1,k);
        M = info(1,2:3,k);
        N = info(2,2:3,k);
        E = info(3,2:3,k);
        S = info(4,2:3,k);
        W = info(5,2:3,k);

        %% Find orthogonal projection of master node
        phi_x = (phi(E(1),E(2))-phi(W(1),W(2)))/(2*hx);
        phi_z = (phi(N(1),N(2))-phi(S(1),S(2)))/(2*hz);
        phi_xx = (phi(E(1),E(2))-2*phi(M(1),M(2))+phi(W(1),W(2)))/hx^2;
        phi_zz = (phi(N(1),N(2))-2*phi(M(1),M(2))+phi(S(1),S(2)))/hz^2;
        phi_xz = (phi(E(1)+1,E(2))-phi(E(1)-1,E(2))-phi(W(1)+1,W(2))+phi(W(1)-1,W(2)))/(4*hx*hz);
        Pnt = [phi_x; phi_z];
        He = [phi_xx phi_xz; phi_xz phi_zz];
        grad = [phi_x, phi_z];
        alpha = -2*phi(node)/(grad*Pnt + sqrt((grad*Pnt).^2 - 4*(Pnt.'*He*Pnt/2)*phi(node)));
        if ~isreal(alpha)
            disp('alpha came back real')
            return
        end
        proj(:,k) = [x(node) + alpha*Pnt(1); z(node) + alpha*Pnt(2)];
        
    end
end