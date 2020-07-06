% https://github.com/augucarv
%
% This code implements the following paper: 
% 
% Li, Y., & Assouar, B. M. (2016). Acoustic metasurface-based perfect 
% absorber with deep subwavelength thickness. Applied Physics Letters, 
% 108(6), 063502.
% _________________________________________________________________________

clear all
close all

%% Geometric parameters [mm] ([mm^2] to S)

d = 3.5e-3;                                                                 % Pore diameter
t = 0.2e-3;                                                                 % Surface thickness
a = 100e-3;                                                                 % UC size
w = 12e-3;                                                                  % Gap between spiral walls
b = 1e-3;                                                                   % Thickness of spiral walls
S = a^2;                                                                    % Cross-sectional area of UC
S_l = w^2;
leff = 635e-3;

%% Frequency domain and fluid parameters

f = 100:0.01:150;
omega = 2*pi.*f;
rho0 = 1.21;                                                                % Air density [kg/m^3]
c0 = 343;                                                                   % Sound speed in air [m/s]                                                                % Velocidade do som no ar [m/s]
z0 = rho0*c0;                                                               % Impedance of air [Rayls]
ni = 1.56e-5;                                                               % Dynamic viscosity [Pa.s]
dv = sqrt(2*ni./(rho0*omega));

%% Theoretical model

K = (d/2).*sqrt(omega./ni);
p = (pi*d^2)/(4*S);
xh = ((32*ni*t./(p*c0*d^2)).*(sqrt(1+(K.^2/32))+(sqrt(2)*K*d)/(8*t)));
yh = ((omega.*t/(p*c0)).*(1+1./sqrt(9+0.5*K.^2)+0.85*(d/t)));
yc = ((-(S/S_l).*cot(omega.*leff/c0)));
alpha = 4*xh./((1+xh).^2 + (yh+yc).^2);

%% Plot

figure()
subplot(3,1,1)
plot(f,alpha,'k','linewidth',1.5)
xlim([min(f) max(f)])
ylim([0 1])
grid on
set(gca,'fontsize',18)
xlabel('Frequency [Hz]')
ylabel('\alpha')
title('Li, Y., & Assouar, B. M. (2016)')
subplot(3,1,2)
plot(f,yh+yc,'b','linewidth',1.5)
xlim([min(f) max(f)])
ylim([-10 10])
grid on
set(gca,'fontsize',18)
xlabel('Frequency [Hz]')
ylabel('ys (yh+yc)')
subplot(3,1,3)
plot(f,xh,'r','linewidth',1.5)
xlim([min(f) max(f)])
ylim([0.5 1.5])
grid on
set(gca,'fontsize',18)
xlabel('Frequency [Hz]')
ylabel('xs (xh)')
