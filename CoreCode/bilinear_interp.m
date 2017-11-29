%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%This Function constructs the bilinear interpolant of the function phi(x,y) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = bilinear_interp(xi,yj,phi,hx,hy)
xbar = @(k,x_i,x) 1+(2*k-1)*((2/hx)*(x-x_i)-1); 
ybar = @(k,y_j,y) 1+(2*k-1)*((2/hy)*(y-y_j)-1);

G = @(x,y) (phi(1)*xbar(0,xi(1),x).*ybar(0,yj(1),y)...
    + phi(2)*xbar(1,xi(1),x).*ybar(0,yj(1),y)...
    + phi(3)*xbar(0,xi(1),x).*ybar(1,yj(1),y)...
    + phi(4)*xbar(1,xi(1),x).*ybar(1,yj(1),y))/4;
end