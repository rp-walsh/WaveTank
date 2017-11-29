% Define solution parameters

%g = 9.81;
g = .001;
%bbeta = 1; % Entropy coeff 
%cp = 2; % above interface
cp = linspace(.1,2,100); % above interface
%cm = 1; % below interface
cm = linspace(.1,2,100); % below interface
H = 1; % height of domain
h = 1/2; % height of the interface
%h = 0.4167; % height of the interface
%h = sqrt(2)/2; % height of the interface
%h = 0.45; % height of the interface
n = 2; % fourier number (even?)
%n = 14; % fourier number (even?)
k = n*pi; % x wave number

[Cp,Cm] = meshgrid(cp,cm);
SolnStructP = zeros(size(Cp));
SolnStructM = zeros(size(Cp));

for ii = 1:length(Cp(:))
    ii
    cp = Cp(ii);
    cm = Cm(ii);
    
    %% Define verticle wave strcture
    alphap = @(k,omega) sqrt(g*cp*k.^2./omega.^2 - k.^2);
    alpham = @(k,omega) sqrt(g*cm*k.^2./omega.^2 - k.^2);
    
    %% Define dispersion relation
    F = @(k,w) alpham(k,w).*sin(alpham(k,w)*h).*cos(alphap(k,w)*(h-H)).*...
        (g*(cp-cm)./(g*cm-w.^2)+1) - alphap(k,w).*sin(alphap(k,w)*(h-H)).*cos(alpham(k,w)*h);
    
    omega = fzero(@(w) F(k,w),sqrt(g/2)-1e-3);
    
    SolnStructP(ii) = g*cp/omega^2 - 1;
    SolnStructM(ii) = g*cm/omega^2 - 1;
    
end

figure(1)
surf(Cp,Cm,SolnStructP,'EdgeColor','none')
xlabel('$c^+$','interpreter','latex','fontsize',18)
ylabel('$c^-$','interpreter','latex','fontsize',18)
title('Solution strcture parameter "+"','interpreter','latex','fontsize',18,'fontweight','bold')

figure(2)
surf(Cp,Cm,SolnStructM,'EdgeColor','none')
xlabel('$c^+$','interpreter','latex','fontsize',18)
ylabel('$c^-$','interpreter','latex','fontsize',18)
title('Solution strcture parameter "-"','interpreter','latex','fontsize',18,'fontweight','bold')