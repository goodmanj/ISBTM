%Data from Jason's File and Goldsby and Kohlsted :Superplastic Defomration
%of Ice
% In Glen-Flow model, Pressures were between 9e5 Pa on the top, and goes to
% 0 on the bottom of the ice model

sigma = .3 ; % MPA
T = 100:1:273; %Kelvins
d = 1  ; % mm

%%%% Nothing needs to be changed here

R = 8.314e-3;            % KJ/(K mol)  Gas constant

p = 1.4; % Grain size exponent?

dis_thresh = 258;% Dislocation Threshold

n_dis = 4.0;% Dislocation exponent
%A_dis_cold = 4.0e5;
A_dis_cold = 1.0e5;
Q_dis_cold = 60;
%A_dis_warm = 6.0e28;
A_dis_warm = 25.0e28;
Q_dis_warm = 180;

gbs_thresh = 255;

n_gbs = 1.8;
A_gbs_cold = 3.9e-3;
Q_gbs_cold = 49;
A_gbs_warm = 3.0e26;
%Q_gbs_warm = 192;
Q_gbs_warm = 190;

n_basal = 2.4;
A_basal = 5.5e7;
Q_basal = 60;
%Diffusion Creep Parameters
b = 4.52e-10 ;% Burger's Vector, in meters
Vm = 1.97e-5; % Molar Volume in m^3
Dov = 9.1e-4; %PreExponential Volume Diffusion in m^2/s
Qv = 59.4; %Grain Boundary Width in kJ/Mol
delta = 9.04e-10; % Grain boundarty Width, Meters
Qb = 49 ; %Activation Energy, Boundary Diffusion in kJ/mol


%%Glen Flow Law:
A0_cold = 3.61e-13;   % Pa^(-3) / s
A0_hot = 1.734e3;     % Pa^(-3) / s
%A0_cold = 7.23e-12;   % Pa^(-3) / s            <--Alternate values
%A0_hot = 3.47e4;     % Pa^(-3) / s
Q_cold = 60e3;        % J/mol
Q_hot = 139e3;        % J/mol
Tcrit = 263.15;       % Use *_cold when T<Tcrit, *_hot when T>Tcrit
R_glen = 8.314;            % J/(K mol)  Gas constant

A = (T<Tcrit) .*A0_cold.*exp(-Q_cold./(R_glen*T)) +  ...
    (T>=Tcrit).*A0_hot.* exp(-Q_hot./(R_glen*T));
% strain rate = A (stress)^n
% stress = 2*mu*strainrate
mu_glen = 1./(2.*A.*sigma^2) ;
%Glen Flow law results in numbers too small for matlab!!!


%%%%%%%%%%%%%%%%%% Dislocation

Z_dis = ((T>dis_thresh).*A_dis_warm.*exp(-Q_dis_warm./(R*T)) + ...
        (T<=dis_thresh).*A_dis_cold.*exp(-Q_dis_cold./(R*T)));
mu_dis = (Z_dis * sigma.^(n_dis-1) ).^-1 ;
    
e_dis = ((T>dis_thresh).*A_dis_warm.*exp(-Q_dis_warm./(R*T)) + ...
        (T<=dis_thresh).*A_dis_cold.*exp(-Q_dis_cold./(R*T))) ...
	.*sigma.^n_dis;
%%%%%%%%%%%%% GBS
Z_gbs = ((T>gbs_thresh).*A_gbs_warm.*exp(-Q_gbs_warm./(R*T)) + ...
        (T<=gbs_thresh).*A_gbs_cold.*exp(-Q_gbs_cold./(R*T)))./d.^p;

mu_gbs = (Z_gbs .* sigma.^(n_gbs-1) ).^-1 ;       
e_gbs = ((T>gbs_thresh).*A_gbs_warm.*exp(-Q_gbs_warm./(R*T)) + ...
        (T<=gbs_thresh).*A_gbs_cold.*exp(-Q_gbs_cold./(R*T))) ...
	.*(sigma.^n_gbs)./d.^p;
%%%%%%% Basal
Z_basal = A_basal.*exp(-Q_basal./(R*T));

mu_basal = (Z_basal .* sigma.^(n_basal -1)).^-1 ;

e_basal = A_basal.*exp(-Q_basal./(R*T)).*(sigma.^n_basal);
d = .001;
% Diffusion Creep
Z_diff = 42 * Vm .*(Dov*exp(-Qv./(R.*T)) + pi*delta * Dov.*exp(-Qb./(R.*T)) /(d))./(R * 1000.*T.*(d)^2); % The 1000 is important, since the GBS wants its units to be in mm, but this wants meters!! In addition, the 'R' is in KJ/K*mol, and needs to be in J/K*mol
mu_diff = (Z_diff .* sigma.^(1 -1)).^-1 ;

%%%other stuff
e = 1./(1./e_basal + 1./e_gbs) + e_dis;

%The plots!
loglog(T,mu_basal,T,mu_dis,T,mu_gbs,T,mu_diff,T,mu_glen)
xlabel('Temperature (K)')
ylabel('Effective Viscosity(Pa s)')
title(['Effective viscosity with Pressure = ' num2str(sigma) 'MPA and Grain size =  ' num2str(d)  'mm'])
legend('T = mu\_basal' , 'T = mu\_dis','T = mu\_gbs' , 'T = mu\_diff','T = Glen Flow Viscosity')
grid