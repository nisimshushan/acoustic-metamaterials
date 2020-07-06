% https://github.com/augucarv
% 
% This code implements the following paper: 
%
% Craggs, A., & Hildebrandt, J. G. (1986). The normal incidence absorption 
% coefficient of a matrix of narrow tubes with constant cross-section. 
% Journal of sound and vibration, 105(1), 101-107.
% _________________________________________________________________________

clear all

%% Frequency
f = 1:1:4000;
omega = 2*pi*f;
%% Material parameters (Table 2, p 6)

P0 = 1.01325e5; % Atmospheric pressure [Pa]
rho0 = 1.23; % Fluid density [kg/m^2]
gama = 1.4; % Specific heat ratio [dim]
Pr = 0.74; % Prandtl Number [dim]
eta = 1.51e-5; % Dynamic viscosity [Pa.s]
nu = eta/rho0; % Kinematic viscous coefficent
nu_ = eta/Pr; % Kinematic thermal coefficent
ce = 343; % Speed of sound in air [m/s]
%% Channel's length
L = 76e-3; % Channel's length
%% Theoretical parameters

r_bar = [1e-3 0.39e-3 0.23e-3 1e-3]; % Hydraulic radius: circ,squar,tri and slit
Beta = sqrt(rho0*omega./eta).*r_bar';
a1 = [8 7 6.5 12];
a2 = [-0.404 -0.513 -0.567 -0.183];
a3 = [0.219 0.218 0.313 0.0938];
a4 = [-8.33e-3 -1.25e-2 -1.46e-2 -1.04e-3];

b1 = [1.33 1.38 1.44 1.20];
b2 = [0.0133 0.0187 0.0252 0.3e-3];
b3 = [-7.5e-3 -1.06e-2 -1.44e-2 -1.88e-3];
b4 = [4.38e-4 6.46e-4 8.96e-4 1.04e-4];

for i = 1:length(a1)

    Re_temp = (eta./r_bar(i)^2).*(a1(i) + a2(i).*Beta(i,:) + a3(i).*Beta(i,:)...
        .^2 + a4(i).*Beta(i,:).^3);
    Re(i) = struct('Re',Re_temp);
    
    rhoe_temp = rho0*(b1(i) + b2(i).*Beta(i,:) + b3(i).*Beta(i,:)...
        .^2 + b4(i).*Beta(i,:).^3);
    rhoe(i) = struct('rhoe',rhoe_temp);
    
    plot1_temp = a1(i) + a2(i).*Beta(i,:) + a3(i).*Beta(i,:).^2 ...
         + a4(i).*Beta(i,:).^3;
    plot1(i) = struct('plot1',plot1_temp);  

end

clear rhoe_temp Re_temp plot1_temp

figure()
plot(Beta(1,:),plot1(1).plot1,Beta(2,:),plot1(2).plot1,Beta(3,:),...
    plot1(3).plot1,Beta(4,:),plot1(4).plot1,'linewidth',2)
xlim([0 10])
legend('Circle','Square','Triangle','Rectangular Slit')
xlabel('\beta')
ylabel({'$R_e \bar{r}^2 / \eta$'},'Interpreter','latex')
set(gca,'FontSize',24)
legend('Location','best')

figure()
plot(Beta(1,:),rhoe(1).rhoe./rho0,Beta(2,:),rhoe(2).rhoe./rho0,...
    Beta(3,:),rhoe(3).rhoe./rho0,Beta(4,:),rhoe(4).rhoe./rho0,'linewidth',2)
xlim([0 10])
legend('Circle','Square','Triangle','Rectangular Slit')
xlabel('\beta')
ylabel({'$\rho_e / \rho_0$'},'Interpreter','latex')
set(gca,'FontSize',24)
legend('Location','best')

figure()
plot(Beta(1,:),Re(1).Re.*(r_bar(1).^2)./eta,Beta(2,:),Re(2).Re.*(r_bar(2).^2)./eta,Beta(3,:),...
    Re(3).Re.*(r_bar(3).^2)./eta,Beta(4,:),Re(4).Re.*(r_bar(4).^2)./eta,'linewidth',2)
xlim([0 10])
legend('Circle','Square','Triangle','Rectangular Slit')
xlabel('\beta')
ylabel({'$R_e \bar{r}^2 / \eta$'},'Interpreter','latex')
set(gca,'FontSize',24)
legend('Location','best')

