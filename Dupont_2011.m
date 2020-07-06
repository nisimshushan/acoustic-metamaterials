% https://github.com/augucarv
%
% This code implements the following paper: 
%
% DUPONT, Thomas et al. Acoustic properties of air-saturated porous 
% materials containing dead-end porosity. Journal of applied physics, 
% v. 110, n. 9, p. 094903, 2011.
% _________________________________________________________________________
clear all
%% Air properties at 20°C

rho0 = 1.205;                                                               % Air density [kg/m^3]
c0 = 343;                                                                   % Sound speed in air [m/s]
P0 = 101325;                                                                % Atmospheric pressure [Pa]
z0 = rho0*c0;                                                               % Characteristic impedance of the air [Rayls]
ni = 15.11e-6;                                                              % Bulk viscosity of air [Pa.s]
Cp = 1.005;                                                                 % Specific heat [kJ/kg.K]
ka = 0.0257;                                                                % Thermal conductivity [W/m.K]
gamma = 1.42;                                                               % Cp/Cv [adimensional]
Pr = 0.713;                                                                 % Prandtl's number
%% Frequency data

f = 1:1:4000;                                                               % Range [Hz]
omega = 2*pi*f;                                                             % Angular frequency [rad/s]
k0 = omega./c0;                                                             % Wave number [1/m]

%% Main pore - Biot's Domain

l = 30e-3;                                                                  % Main pore's length [cm]
d = 2e-3;                                                                   % Main pore's radius [mm]
dmin = 1.8e-3;
phiB = 14/100;
sigmaB = (32*ni)/(phiB*d^2);                                                % Flux resistivity [Pa.s/N]
LambdaB = dmin/2;                                                           % Viscous characteristic length [mm]
LambdalB = d/2;                                                             % Thermal characteristic length [mm]
omegacB = sigmaB*phiB/rho0;                                                 % Eq. (16)

FB = sqrt(1+1i.*((omega.*4*ni*rho0*1)./((phiB^2)*(sigmaB^2)*(LambdaB^2)))); % Eq. (17)
GB = sqrt(1+1i.*((omega.*rho0*Pr*LambdalB^2)./(16*ni)));                    % Eq. (18)
rhoB = (rho0/phiB)*(1-1i*(omegacB./omega).*FB);                             % Eq. (14)
KB = (1/phiB)*((gamma*P0)./(gamma-(gamma-1)*...
    (1-1i*((8*ni.*GB)./(omega.*Pr*rho0*LambdalB^2))).^(-1)));               % Compressibility - Biot - Eq. (15)
ZB = sqrt(rhoB.*KB);                                                        % Impedance - Biot - Eq. (19)
kB = omega.*sqrt(rhoB./KB);                                                 % Wavenumber - Biot - Eq. (19)

%% DE pores - Dupont's Domain

lDE = 25e-3;                                                                % DE pore's length [cm]
phiDE = 13.5/100;
sigmaDE = (32*ni)/(phiDE*d^2);                                              % Flux resistivity [Pa.s/N]
omegacDE = sigmaDE*phiDE/rho0;                                              % Eq. (16)

FDE = sqrt(1+1i*((omega.*4*ni*rho0*1)./((phiDE^2)*(sigmaDE^2)*...
    (LambdaB^2))));                                                         % F - Eq. (17)
GDE = sqrt(1+1i.*((omega.*rho0*Pr*LambdalB^2)./(16*ni)));                   % G - Eq. (18)
rhoDE = (rho0/phiDE)*(1-1i*(omegacDE./omega).*FDE);                         % Eq. (14)
KDE = (1/phiDE)*((gamma*P0)./(gamma-(gamma-1)*...
    (1-1i*((8*ni*GDE)./(omega.*Pr*rho0*LambdalB^2))).^(-1)));               % Compressibility - Eq. (15)
ZDE = sqrt(rhoDE.*KDE);                                                     % Impedance - Eq. (19)
kDE = omega.*sqrt(rhoDE./KDE);                                              % Wavenumber - DE - Eq. (19)

%% Transfer matrices

for i = 1:length(f)
   TDE_temp = [cos(kDE(i).*lDE)    1i*ZDE(i).*sin(kDE(i).*lDE)
               (1i./ZDE(i))*sin(kDE(i).*lDE)   cos(kDE(i)*lDE)];
   TDE(i) = struct('TDE', TDE_temp);                                        % Transfer matrix DE - Eq. (21)
end

for i = 1:length(f)
   YDE_temp = [(TDE(i).TDE(2,2)./TDE(i).TDE(1,2))    -1./TDE(i).TDE(1,2)
               -1./TDE(i).TDE(1,2)     (TDE(i).TDE(1,1)./TDE(i).TDE(1,2))];   
   YDE(i) = struct('Y', YDE_temp);                                          % Transfer matrix DE - Eq. (23)
end

for i = 1:length(f)
   TB_temp = [cos(kB(i).*l)    1i*ZB(i).*sin(kB(i).*l)
               (1i./ZB(i))*sin(kB(i).*l)   cos(kB(i)*l)];
   TB(i) = struct('TB', TB_temp);                                           % Transfer matrix Biot - Eq. (26)
end

for i = 1:length(f)
   YB_temp= [(TB(i).TB(2,2)./TB(i).TB(1,2))    -1./TB(i).TB(1,2)
               -1./TB(i).TB(1,2)     (TB(i).TB(1,1)./TB(i).TB(1,2))];   
   YB(i) = struct('Y', YB_temp);                                            % Transfer matrix Biot - Eq. (28)
end

for i = 1:length(f)
    
    TNS_temp = (1./YB(i).Y(1,2)).*...
               [-YB(i).Y(2,2)      -1
               ((YB(i).Y(1,2)).^2)-(YB(i).Y(2,2).*((YB(i).Y(1,1))+...
               (YDE(i).Y(1,1))-(((YDE(i).Y(1,2).^2)./...
               ((YDE(i).Y(2,2)))))))   -((YB(i).Y(1,1))+...
               (YDE(i).Y(1,1))-(((YDE(i).Y(1,2).^2)./...
               ((YDE(i).Y(2,2))))))];
    
    TNS(i) = struct('TNS', TNS_temp);                                       % Transfer matrix - Eq. (31)
    
end

clearvars TDE_temp YDE_temp TB_temp YB_temp TNS_temp

%% Absorption coefficient - Material with hardback at the end

for i = 1:length(f)
ZsHB_temp = TNS(i).TNS(1,1)./TNS(i).TNS(2,1);                               % Eq. (45)
ZsHB(i) = struct('ZsHB',ZsHB_temp);
alphaHB(i) = 1-abs((ZsHB(i).ZsHB-z0)./(ZsHB(i).ZsHB+z0)).^2;
end

%% Absorption coefficient  - Cavity - Eq. (44)

lcav = 20e-3;

for i = 1:length(f)
    TCAV_temp = [cos(k0(i).*lcav)  1i*z0.*sin(k0(i)*lcav)
                 (1i./z0)*(sin(k0(i)*lcav))  cos(k0(i)*lcav)];
    TCAV(i) = struct('TCAV', TCAV_temp);
end

for i = 1:length(f)
    TS_temp = TNS(i).TNS*TCAV(i).TCAV;
    TS(i) = struct('TS', TS_temp);
end

for i = 1:length(f)
ZsCAV_temp = TS(i).TS(1,1)./TS(i).TS(2,1);                                  % Eq. (45)
ZsCAV(i) = struct('ZsCAV',ZsCAV_temp);
alphaCAV(i) = 1-abs((ZsCAV(i).ZsCAV-z0)./(ZsCAV(i).ZsCAV+z0)).^2;
end

plot(f,alphaCAV,f,alphaHB,'k','linewidth',1.5)
set(gca,'fontsize',14)
xlabel('Frequ�ncia, kHz')
ylabel('Coeficiende de absor��o \alpha')
grid on