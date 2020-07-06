% https://github.com/augucarv
%
% This code implements the following paper:
%
% SALISSOU, Yacoubou; PANNETON, Raymond; DOUTRES, Olivier. Complement to 
% standard method for measuring normal incidence sound transmission loss with 
% three microphones. The Journal of the Acoustical Society of America, 
% v. 131, n. 3, p. EL216-EL222, 2012.
% 
% The virtual measurement was made in COMSOL software.

% _________________________________________________________________________
clear all
%% Air properties at 20°C
f = 1:1:1000;                                                               % Frequency [Hz]
omega = 2*pi*f;                                                             % Angular frequency [rad/s]
c0 = 343;                                                                   % Speed of sound [m/s]
rho0 = 1.205;                                                               % Air density [kg/m^3]
k = omega/c0;                                                               % Wavenumber [1/m]
Z0 = rho0*c0;                                                               % Characteristic impedance [Rayls]

%% Tube's dimensions

s =  18e-3;                                                                 % Distance between Mic 1 and Mic 2
Da = 25e-3;                                                                 % Air cavity's length
l = 22e-3;                                                                  % Distance between Mic 2 and sample
d = 20e-3;                                                                  % Sample's length
L  =d+Da;

%% Transfer functions
load('Microfones.mat')
H12 = Microfones(:,2)./Microfones(:,1);                                 % Transfer function Mic2/Mic1
H13 = Microfones(:,3)./Microfones(:,1);                                 % Transfer function Mic3/Mic1

%% Transfer matrix

pa_0 = -2*1i*exp(1i*k*l).*(H12'.*sin(k*(l+s))-sin(k*l))./...
    (H12'.*exp(-1i*k*s)-1);
ua_0 = ((2*exp(1i*k*l))/Z0).*(H12'.*cos(k*(l+s))-cos(k*l))./...
    (H12'.*exp(-1i*k*s)-1);
pa_d = -2*1i*exp(1i*k*l).*(H13'.*sin(k*s).*cos(k*Da))./...
    (H12'.*exp(-1i*k*s)-1);
ua_d = ((2*exp(1i*k*l))/Z0).*(H13'.*sin(k*s).*sin(k*Da))./...
    (H12'.*exp(-1i*k*s)-1);

T = [(pa_d.*ua_d+pa_0.*ua_0)./(pa_0.*ua_d+ua_0.*pa_d) ;                     % T11
     (pa_0.^2-pa_d.^2)./(pa_0.*ua_d+ua_0.*pa_d) ;                           % T12
     (ua_0.^2 -ua_d.^2)./(pa_0.*ua_d+ua_0.*pa_d) ;                          % T21
     (pa_d.*ua_d+pa_0.*ua_0)./(pa_0.*ua_d+ua_0.*pa_d)];                     % T22

%% Plotting
R = (T(1,:)-Z0*T(3,:))./(T(1,:)+Z0*T(3,:));
alpha = 1-abs(R.^2);
plot(f,alpha)