% https://github.com/augucarv
% 
% This code implements the following paper: 
%
% Huang, S., Fang, X., Wang, X., Assouar, B., Cheng, Q., & Li, Y. (2018). 
% Acoustic perfect absorbers via spiral metasurfaces with embedded apertures.
% Applied Physics Letters, 113(23), 233501.
% _________________________________________________________________________
clear all

xopt = [0.0088    0.0074    0.0386];

%% Frequency
f = 1:0.1:1000;
omega = 2*pi*f; 
%% Frequency domain and fluid parameters

Temp = 293.15;                                                              % Temperatura
mu = -8.38278E-7+8.35717342E-8*Temp^1-7.69429583E-11*Temp^2+...
    4.6437266E-14*Temp^3-1.06585607E-17*Temp^4;                             % Dynamic viscosity of air
Cp = 1047.63657-0.372589265*Temp^1+9.45304214E-4*Temp^2-...
    6.02409443E-7*Temp^3+1.2858961E-10*Temp^4;                              % Calor espec�fico [kJ/kg.K]
K = -0.00227583562+1.15480022E-4*Temp^1-7.90252856E-8*Temp^2+...
    4.11702505E-11*Temp^3-7.43864331E-15*Temp^4;                            % Condutividade t�rmica
P0 = 101325;  

rho0 = 1.21;                                                                % Air density [kg/m^3]
c0 = sqrt(1.4*287*(Temp));                                                  % Sound speed in air [m/s]                                                                % Velocidade do som no ar [m/s]
k0 = omega./c0;                                                             % Wave Number
Z0 = rho0*c0;                                                               % Impedance of air [Rayls]
gamma = 1.4;

%% Geometry parameters

h = xopt(1); % Cross section's height
w = xopt(2); % Cross section's width
af = xopt(3);% Outer radius of spiral

a = 40e-3; % UC's side length
A = a^2; % UC's total area

s0 = w*h; % Cross section area of the entrance of the channel

gap = 1e-3;                                                                 % Wall thickness
ai = w/2;                                                                   % Starting radius
b = (w+gap)/(2*pi);                                                         % Spiral growth ratio
n = (af-ai)/(2*pi*b);                                                       % Number of turns
theta0 = 0;                                                                 % Initial theta
theta_final = (af-ai)/b + theta0;                                           % Final theta
syms theta
fun = sqrt((ai + b*theta).^2 + b^2);
L = double(int(fun,[theta0 theta_final]));                             % Length of spiral

%% Theoretical parameters

kv = sqrt(-1j*omega*rho0/mu); % Viscous wave number (2)
kh = sqrt(-1j*omega*rho0*Cp/K); % Thermal wavenumber (3)

Phi_v = 0;
Phi_h = 0;

for m = 0:100

m_ = (m+1/2)*pi;
alpha_mv = sqrt(kv.^2-(2*m_/w).^2);
beta_mv = sqrt(kv.^2-(2*m_/h).^2);
alpha_mh = sqrt(kh.^2-(2*m_/w).^2);
beta_mh = sqrt(kh.^2-(2*m_/h).^2);

Phi_vtemp = (kv.^2).*((((alpha_mv*m_).^(-2)).*(1-...
    (tan(alpha_mv*w/2)./(alpha_mv*w/2)))+...
    ((((beta_mv*m_).^(-2)).*(1-...
    (tan(beta_mv*h/2)./(beta_mv*h/2))))))); % Function of viscous field (A15)

Phi_htemp = (kh.^2).*((((alpha_mh*m_).^(-2)).*(1-...
    (tan(alpha_mh*w/2)./(alpha_mh*w/2)))+...
    ((((beta_mh*m_).^(-2)).*(1-...
    (tan(beta_mh*h/2)./(beta_mh*h/2))))))); % Function of thermal field (A16)

    Phi_v = Phi_v + Phi_vtemp;
    Phi_h = Phi_h + Phi_htemp;
    
 end

rhoc = rho0./Phi_v; % Complex air density (7)
cc = c0*sqrt(Phi_v./(gamma*(gamma-1).*Phi_h)); % Complex velocity (COMSOL Equation)
kc = k0.*sqrt((gamma-(gamma-1)*Phi_h)./Phi_v); % Complex wavenumber (6)
Z = -1j*rhoc.*cc.*(A./(s0)).*cot(kc*L); % Impedance of the system (A20)
R = abs((Z-Z0)./(Z+Z0)); % Reflection coefficient (A19)
alpha = 1-abs(R).^2; % Sound absorption coefficient (1)

figure()
% subplot(2,1,1)
plot(f,alpha,f,abs(R),'--k','linewidth',2)
xlim([min(f) max(f)])
ylim([0 1])
grid on
set(gca,'fontsize',24)
xlabel('Frequency [Hz]')
ylabel('|R|^2')

    