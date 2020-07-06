% https://github.com/augucarv
%
% This code implements the following paper: 
%
% LECLAIRE, Philippe et al. Acoustical properties of air-saturated porous 
% material with periodically distributed dead-end pores a. The Journal of 
% the Acoustical Society of America, v. 137, n. 4, p. 1772-1782, 2015.
% _________________________________________________________________________
clear all
f = 100:1:1000;                                                             % Frequency range [Hz]

%% Air properties at 20°C

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
z0 = rho0*c0;                                                               % Characteristic impedance of air [Rayls]
ni = 1.511e-5;                                                              % Bulk viscosity [Pa.s]
gamma = 1.42;                                                               % Cp/Cv [adimensional]
Pr = (ni*Cp)/ka;                                                            % Prandtl's number

%% Dimensions

amp = 2e-3;                                                                 % Main pore's radius [mm]
ade = 1.6e-3;                                                               % DE pore's radius [mm]
Lmp = 4.5e-2;                                                               % Main pore's length [cm]
d = 0.85e-2;                                                                % DE pore's length [cm]

h = 4.5e-3;                                                                 % Spacing between DE's pores [cm]
N = 8;                                                                      % Number of DE pores per node
n = 10;                                                                     % Number of full periods

omega = 2*pi*f;                                                             % Angular frequency [rad/s]

%% Caracter�sticas do poro principal (Biot's Domain)

Amp = pi*amp^2;                                                             % Main pore's cross section area [mm^2]
sigma_mp = (8*ni)/(amp^2);                                                  % Flux resistivity [Pa.s/N]
ka_mp = (amp^2)/8;                                                          % Thermal permeability [W/mK]
Lambda_mp = amp;                                                            % Viscous characteristic length [mm]
Lambdal_mp = amp;                                                           % Thermal characteristic length [mm]
omegal_mp = omega.*sqrt(Pr);
omegab_mp = ((sigma_mp^2)*(Lambda_mp^2))./(4*1*rho0*ni);
omegalb_mp = (Lambdal_mp^2)./(4*rho0*ka_mp^2);
rho_mp = rho0*1*(1+(sigma_mp./(-1i*omega*rho0)).*...                        % Effective density 
    sqrt(1+(-1i*omega./omegab_mp)));                                        % Eq. (38)
C_mp = (1/(rho0*c0^2))*(gamma-((gamma-1)./...                               % Effective compressibility
    (1+(ni./(-1i*omegal_mp*rho0*ka_mp)).*...                                % Eq. (39)
    sqrt(1+(-1i*omegal_mp/omegalb_mp)))));
Zmp = sqrt(rho_mp./C_mp);                                                   % Impedance
kmp = omega.*sqrt(rho_mp.*C_mp);                                            % Dissipative wavenumber

%% DE pores

Ade = pi*ade^2;                                                             % Cross section area [mm^2]
sigma_de = 8*ni/(ade^2);                                                    % Flux resistivity [Pa.s/N]
ka_de = (ade^2)/8;                                                          % Thermal permeability [W/mK]
Lambda_de = ade;                                                            % Viscous characteristic length [mm]
Lambdal_de = ade;                                                           % Thermal characteristic length [mm]
omegal_de = omega*sqrt(Pr);
omegab_de = (sigma_de^2)*(Lambda_de^2)/(4*1*rho0*ni);
omegalb_de = (Lambdal_de^2)/(4*rho0*ka_de^2);
rho_de = rho0*1*(1+(sigma_de./(-1i*omega*rho0)).*...                        % Effective density
    sqrt(1+(-1i*omega./omegab_de)));
C_de = (1/(rho0*c0^2))*(gamma-((gamma-1)./...                               % Effective compressibility
    (1+(ni./(-1i*omegal_de*rho0*ka_de)).*...
    sqrt(1+(-1i*omegal_de/omegalb_de)))));
Zde = sqrt(rho_de./C_de);                                                   % Effective impedance
kde = omega.*sqrt(rho_de.*C_de);                                            % Dissipative wavenumber

%% TMM

X = 1i*(N/2)*(Ade/Amp)*(Zmp/Zde).*tan(kde.*d);                              % Periodicity impedance - Eq. (4)
qh = (acos(cos(kmp*h)+1i*X.*sin(kmp*h)))./h;                                % Bloch's wavenumber - Eq. (1)

%% LFA

q = kmp.*sqrt(1+((2*X)./(1i*kmp*h)));                                       % Bloch's wavenumber - Eq. (32)

%% Sound absorption coefficient - TMM

y = exp(1i*kmp*h);                                                          % Eq. (6)

for i = 1:length(f)
Tc_temp = [(1+X(i))*y(i)  X(i)                                              % Eq. (7)
          -X(i)  (1-X(i))/y(i)];
Tc(i) = struct('Tc',Tc_temp);
M(i) = struct('M',Tc_temp^n);                                               % Eq. (8)
end

for i = 1:length(f)
rn_inf(i) = -(M(i).M(2,1))./(M(i).M(2,2));                                  % reflection coefficient, infinite main pore - Eq. (10a)
alpha_inf(i) = 1-abs(rn_inf(i)).^2;                                         % absorption coefficient, infinite main pore - Eq. (19)
rn_hb(i) = (M(i).M(1,1)-M(i).M(2,1))/(M(i).M(2,2)-M(i).M(1,2));             % reflection coefficient, hardbacked main pore - Eq. (12)
alpha_hb_TMM(i) = 1-abs(rn_hb(i)).^2;                                       % absorption coefficient, hardbacked main pore - Eq. (19)
end

clearvars  Tc_temp
%% SLAB

A = pi*(13.25e-3)^2;

phi = Amp/A;
                                                                            % Perforation rate - Eq. (13)
phil = phi*(z0./Zmp);                                                           

for i = 1:length(f)
   T_temp = [((1+phil(i))/(2*phi))       -((1-phil(i))/(2*phi))             % Eq. (15)
            -((1-phil(i))/(2*phil(i)))   ((1+phil(i))/(2*phil(i)))];
   T(i) = struct('T',T_temp);
end

for i = 1:length(f)
   Mil_temp = M(i).M*T(i).T;                                                % Eq. (18)
   Mil(i) = struct('Mil',Mil_temp);
end

for i =1:length(f)
    Rn(i) = -Mil(i).Mil(2,1)./Mil(i).Mil(2,2);                              % ref. coef. open slab - Eq. (16a)
    Rnil(i) = (Mil(i).Mil(1,1)-Mil(i).Mil(2,1))./...                        % ref. coef. hb slab - Eq. (17)
        (Mil(i).Mil(2,2)-Mil(i).Mil(1,2));
    alpha_hb_slab(i) = 1-abs(Rnil(i).^2);                                   % abs. coef. hardbacked main pore - Eq. (19)
end
clearvars  T_temp Mil_temp

%% Plot

figure()

plot(f,alpha_hb_slab,'-k','linewidth',2)
set(gca,'fontsize',18)
xlabel('Frequency [Hz]')
ylabel('Sound absorption coefficient, \alpha')
grid on
xlim([min(f) max(f)])
ylim([0 1])
title(['L = ' num2str(Lmp*100) ' cm, N = ' num2str(N) ', d = ' num2str(d*100) ' cm, a_{mp} = ' num2str(amp*1000) ' mm, a_{de} = ' num2str(ade*1000) ' mm, ' 'h = ' num2str(h*1000) ' mm, \phi_s = ' num2str(phi*100) ' %'])