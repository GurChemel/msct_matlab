clear; clc;

E_ev = 50e3; % ev
theta = linspace(0,pi,1800);

% 1 ev = 1.602e-12 cm^2 g /s =  1.602e-15 cm^2 kg /s 
% 1 erg = 1
h = 6.626e-34; %m^2kg /s
m = 9.11e-31; %kg
c = 29979245800; %cm/s
r = 2.8e-13; %cm
% m*c*c = cm^2 kg / s^2
E_arg = E_ev*(1.602e-15);

epsilon = 1./(1+E_arg*(1-cos(theta))/(m*c*c));
epsilon_inv = (1+E_arg*(1-cos(theta))/(m*c*c));

Z = 8;
ds_de = pi*r*(m*c*c/E_arg)*Z.*(epsilon_inv+epsilon).*(1-epsilon.*((sin(theta)).^2)./(1+epsilon.^2));
de_dt = -(E_arg/(m*c*c)).*sin(theta)./((1+(E_arg/(m*c*c))*cos(theta)).^2);

ds_dt = abs(ds_de.*de_dt);
polarplot(theta,ds_dt);

E_out = E_arg./(1+(E_arg/(m*c*c))*(1-cos(theta)));
E_out_ev = E_out/(1.602e-15);

%%
clear; clc;
theta = linspace(0,2*pi,1800);

m = 9.11e-31; %kg
c = 29979245800; %cm/s
r = 2.8e-13; %cm

E0 = (50e3)*(1.602e-15);

epsilon = E0/(m*c*c);
sigma = ((1/(2*(r^2)))*(1+(cos(theta)).^2)).*(1+(4*epsilon*epsilon*(sin(theta/2)).^4)./((1+(cos(theta)).^2).*(1+2*epsilon*(sin(theta/2)).^2)))./((1+2*epsilon*(sin(theta/2)).^2).^2);

polarplot(sigma)