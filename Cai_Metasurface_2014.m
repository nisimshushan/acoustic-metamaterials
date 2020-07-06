% https://github.com/augucarv
% 
% This code implements the paper:
% 
% Cai, X., Guo, Q., Hu, G., & Yang, J. (2014). Ultrathin low-frequency sound 
% absorbing panels based on coplanar spiral tubes or coplanar Helmholtz 
% resonators. Applied Physics Letters, 105(12), 121901.
% _________________________________________________________________________

clear all
close all

%% Frequency parameters

f = 200:1:600;
omega = 2*pi.*f;

%% Geometric parameters [mm] ([mm^2] to S)

Ro = 30e-3;                                                                 % Outer radius
rt = 4.85e-3;                                                               % Optimal tube radius
thick = 2*rt;                                                             % Optimal thickness of air channels
xi = (rt^2)/(Ro^2);                                                         % Porosity of panel
gap = 3.4e-3;                                                              % Wall thickness
ai = thick/2;                                                               % Starting radius
b = (thick+gap)/(2*pi);                                                     % Spiral growth ratio
n = (Ro-ai)/(2*pi*b);                                                       % Number of turns
theta0 = 0;                                                                 % Initial theta
theta_final = (Ro-ai)/b + theta0;                                           % Final theta
syms theta
fun = sqrt((ai + b*theta).^2 + b^2);
Lt = double(int(fun,[theta0 theta_final]));
clear theta

%% Air properties at 20C

Temp = 273.15;                                                              % Temperature [K]
eta = -8.38278E-7+8.35717342E-8*Temp^1-7.69429583E-11*Temp^2+...
    4.6437266E-14*Temp^3-1.06585607E-17*Temp^4;
Cp = 1047.63657-0.372589265*Temp^1+9.45304214E-4*Temp^2-...
    6.02409443E-7*Temp^3+1.2858961E-10*Temp^4;                              % Specific heat [kJ/kg.K]
ka = -0.00227583562+1.15480022E-4*Temp^1-7.90252856E-8*Temp^2+...
    4.11702505E-11*Temp^3-7.43864331E-15*Temp^4;                            % Thermal conductivity
P0 = 101325;                                                                % Atmospheric pressure [Pa]

rho0 = 1.204;                                                               % Air density [kg/m^3]
c0 = sqrt(1.4*287*(Temp));                                                  % Speed of sound [m/s]                                                                % Velocidade do som no ar [m/s]
z0 = rho0*c0;                                                               % Characteristic impedance [Rayls]
ni = 1.511e-5;                                                              % Bulk viscosity [Pa.s]
Pr = (ni*Cp)/ka;                                                            % Prandtl's number

Z = rho0*c0;                                                                % Impedance of air [Rayls]
nil = 0.3e-4;
gamma = 1.4;

Carg = rt*(-1i*omega/ni).^(1/2);
rhoarg = rt*(-1i*omega*gamma/nil).^(1/2);
G_rho = besselj(1,rhoarg)./besselj(0,rhoarg); 
G_C = besselj(1,Carg)./besselj(0,Carg); 
%% Theoretical model
                                         
rho = rho0*(1-2*((-1i*omega/ni).^(-1/2)).*G_rho/rt).^(-1);
C = (1/(gamma*P0))*(1+(2*(gamma-1)*(-1i*omega*gamma/nil).^(-1/2)).*...
    G_C/rt);
k = omega.*sqrt(rho.*C);
Zc = sqrt(rho./C);
Zt = -1i*Zc.*cot(k*Lt);
Zin = Zt/xi;
R = (Zin-Z)./(Zin+Z);
alpha = 1-abs(R).^2;

%% Plot

figure()
plot(f,alpha,'b','linewidth',3)
xlim([min(f) max(f)])
ylim([0 1])
grid on
set(gca,'fontsize',28)
xlabel('Frequency [Hz]')
ylabel('\alpha')
%title('Cai et al. (2014)')

