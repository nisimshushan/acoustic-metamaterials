% https://github.com/augucarv
%
% This code implements the following paper: 
%
% S Y Song et al 2017 J. Phys. D: Appl. Phys. 50 235303
% _________________________________________________________________________

clear all

%% Frequency
f = 1:1:4000;
omega = 2*pi*f;
%% Material parameters (Table 2, p 6)

P0 = 1.01325e5; % Atmospheric pressure [Pa]
rho0 = 1.23; % Fluid density [kg/m^2]
gamma = 1.4; % Specific heat ratio [dim]
Pr = 0.74; % Prandtl Number [dim]
mu = 1.79e-5; % Dynamic viscosity [Pa.s]
nu = mu/rho0; % Kinematic viscous coefficent
nu_ = nu/Pr; % Kinematic thermal coefficent

%% Roughness and model properties

a = 1e-3; % Diameter of slit
d = 50e-3; % Length of slit
epsilon = 0; % Relative roughness
beta = 2*pi; % Wave number of roughness

%% Equivalent fluid parameters

sigma0 = 12*mu/a^2; % Flow resistivity in a smooth slit
sigma = sigma0*((((1+2*epsilon^2)./(1-4*epsilon^2)^(2.5)) - ...
    1/(1-2*epsilon)^3).*(2*exp(-beta/(5*pi))/(1+exp(-beta/(5*pi))))+ ...
    1/(1-2*epsilon)^3); % Flow resistivity in rough slit (Eq. 30)
alpha_inf = 1+((exp(beta) + 1)/(exp(beta) - 1)).*beta*epsilon^2; % Tortuosity (Eq. 46)
%% Acoustic parameters (Table 1)

Lambda = ((8*mu*alpha_inf./sigma).^(1/2))./sqrt(2/3); % Characteristic viscous length
func = @(t) sqrt(1+(epsilon.^2).*(beta.^2).*sin(t).^2);
A = integral(func,0,pi/2);
Lambda_ = a*pi./(2.*A); % Characteristic thermal length 
k0 = mu./sigma; % Static viscous permeability
k0_ = (Lambda_^2)/12; % Static thermal permeability

%% Effective density and compressibility - Approximed solutions (Eqs. 47 and 48)
c = 5/6;
rho = (alpha_inf + (nu./(1j.*omega.*k0)).*((1-c+(c.*sqrt((1+...
    ((1j.*omega.*(2.*alpha_inf.*k0./(c*Lambda)).^2)./nu))))))).*rho0; % Effective density (Eq. 47)
C = (gamma - (gamma-1)./(1+(nu_./(1j.*omega.*k0_)).*(1+...
    (1j.*omega./nu_).*(Lambda_/6).^2).^(1/2)))./(gamma*P0); % Effective compressibility (Eq. 48)


%% Effective density and compressibility - Exact solutions (Eqs. 52 and 53)
s = sqrt((omega.*rho0*a^2)./(4*mu));
rho_e = (1./(1-(tanh(s.*sqrt(1j))./(s.*sqrt(1j))))).*rho0; % Effective density (Eq. 52)
C_e = (1+(gamma-1).*(tanh(sqrt(Pr).*s.*sqrt(1j))...
    ./(sqrt(Pr).*s.*sqrt(1j))))./(gamma*P0); % Effective compressibility (Eq. 53)


%% Plot (to reproduce Fig. 7 of the paper)

figure()
subplot(2,2,1)
semilogx((f.*rho0*a^2)./(4*mu),real(rho_e./rho0),(f.*rho0*a^2)./(4*mu),real(rho./rho0),...
        '--k','linewidth',1.5)
xlabel('f\rho_0a^2/4\mu')
ylabel('Real(\rho_0/\rho)')
set(gca,'FontSize', 18)
legend('Exact Solution','Paper')
subplot(2,2,2)
semilogx((f.*rho0*a^2)./(4*mu),imag(rho_e./rho0),(f.*rho0*a^2)./(4*mu),imag(rho./rho0),...
    '--k','linewidth',1.5)
xlabel('f\rho_0a^2/4\mu')
ylabel('Imag(\rho_0/\rho)')
set(gca,'FontSize', 18)
legend('Exact Solution','Paper')
subplot(2,2,3)
semilogx((f.*rho0*a^2)./(4*mu),real(C_e*gamma*P0),(f.*rho0*a^2)./(4*mu),...
    real(C*gamma*P0),'--k','linewidth',1.5)
xlabel('f\rho_0a^2/4\mu')
ylabel('Real(C\times\gammaP_0)')
set(gca,'FontSize', 18)
legend('Exact Solution','Paper')
subplot(2,2,4)
semilogx((f.*rho0*a^2)./(4*mu),imag(C_e*gamma*P0),(f.*rho0*a^2)./(4*mu),...
    imag(C*gamma*P0),'--k','linewidth',1.5)
xlabel('f\rho_0a^2/4\mu')
ylabel('Imag(C\times\gammaP_0)')
set(gca,'FontSize', 18)
legend('Exact Solution','Paper')