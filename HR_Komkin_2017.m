% https://github.com/augucarv
%
% This code implements the following paper: 
%
% Komkin, A. I., M. A. Mironov, and A. I. Bykov. "Sound absorption by a 
% Helmholtz resonator." Acoustical Physics 63.4 (2017): 385-392.
% _________________________________________________________________________
clear all
f = 1:0.1:1000;

%% HR %%
%% Geometrical parameters
% Neck
l = 10.6e-3; %Length
d0 = 3.7e-3; %Diameter
S0 = pi*(d0/2)^2; %Cross section area
% Cavity
L = 13.1e-3; %Length
D = 37.8e-3; %Diameter
S = pi*(D/2)^2; %Cross section area
%V = S*L; %Volume
V = (4/3)*pi*(D/2)^3;
% Derived parameters
Ss = pi*D*L+((pi*D^2)/2)-pi*(d0^2)/4; %Surface area
g = d0/D;
m = S/S0;
%% Frequency viscothermal properties
c = 343; %Speed of sound in air
rho0 = 1.21; % Air density
Z0 = rho0*c; % Air impedance
omega = 2*pi.*f; %Angular frequency
k = omega./c; %Wavenumber
nu = 1.5e-5; %Kinematic viscosity of air
deltav = (2*nu./omega).^(1/2); %Depth of viscous boundary layer
X = 2.1e-5; %Thermal diffusity of medium
deltaX = (2*X./omega).^(1/2); %Depth of thermal boundary layer
%% Theory's parameters 1
la = 0.82*(1-1.34*g)*d0; %Attached length
Deltav = 2*deltav./d0;
le = l.*(1+Deltav)+la; %Effective length of ressonator's neck
omega0 = c*g./(sqrt(L*(le+(L*g^2)/3))); %Corrected natural frequency
k0 = omega0./c; %Corrected wavenumber
%% Theory's parameters 2
N = 0.28;
E = 0.85;
Rv = 2*k.*deltav*(((l/d0)+N+E)); %Normalized viscous resistance
Rx = (S0*Ss.*deltaX./(k.*V^2)); %Normalized thermal resistance
Omegak = (omega./omega0); %Corrected dimensionless frequency
Z_res = (m.*(Rv+Rx)+1i.*(S./(V.*k0)).*(Omegak-(1./Omegak))); %Impedance of resonator
Z_hrr = Z0.*Z_res; %Surface impedance
R = (Z_hrr-Z0)./(Z_hrr+Z0);
alpha = 1-abs(R).^2;

% figure()
% subplot(2,1,1)
% plot(f,abs(R),'k','linewidth',2)
% xlim([min(f) max(f)])
% ylim([0 1])
% grid on
% set(gca,'fontsize',24)
% xlabel('Frequency [Hz]')
% ylabel('|R|^2')
% subplot(2,1,2)
figure()
plot(f,alpha,'b','linewidth',2)
xlim([min(f) max(f)])
ylim([0 1])
grid on
set(gca,'fontsize',24)
xlabel('Frequency [Hz]')
ylabel('\alpha')