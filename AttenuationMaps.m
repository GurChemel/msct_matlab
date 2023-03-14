clear; clc;
addpath('PhotonAttenuation\')

figure(1); clf;
E = logspace(log10(0.001), log10(20), 500);  % define energy grid
mac  = PhotonAttenuationQ(92, E, 'mac');
meac = PhotonAttenuationQ(92, E, 'meac');
loglog(E, mac); hold on;
loglog(E, meac, 'b-.');
legend({'mac', 'meac'});
ylabel('Attenuation in cm^2/g');
xlabel('Photon Energy in MeV');
title({'Photon Attenuation Coefficients for Uranium',...
'see http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z92.html'});

%%
figure(1); clf;
Z = 1:92; % elements with Z in 1-92 range (Elements higher than 92 are not defined)
% Z = [1,6,7,8,15,20];
E = logspace(log10(0.001), log10(0.150), 500);  % define energy grid
meac = PhotonAttenuationQ(Z, E, 'meac');
imagesc(log10(meac)); colorbar;
colormap('jet')
title('Log of Photon Mass Energy-Absorption Coefficients in cm^2/g');
xlabel('Atomic Number of Elements');
ylabel('Energy in KeV');
set(gca,'YTick',linspace(1, length(E), 10));
set(gca,'YTickLabel',round(1e3*logspace(log10(0.001), log10(0.150), 10)))
% set(gca,'XTick',linspace(1, length(E), 10));
set(gca,'XTickLabel',Z)

%%
figure(1); clf;
Z = [1,6,7,8,15,20];
E = logspace(log10(0.001), log10(0.150), 500);  % define energy grid
mac = PhotonAttenuationQ(Z, E, 'mac');
loglog(E, mac); hold on;
legend({'H','C','N','O','P','Ca'});
ylabel('Attenuation in cm^2/g');
xlabel('Photon Energy in MeV');
title({'Photon Attenuation Coefficients for Uranium',...
'see http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z92.html'});

%%

Z = {'COMPACT BONE','Blood','Soft Tissue','Water'};
E = logspace(log10(0.001), log10(0.150), 500);  % define energy grid
mac = PhotonAttenuation(Z, E, 'mac');

P = PhysProps(Z);

mac_no_water   = mac(:,1:(end-1));
mac_only_water = repmat(mac(:,end),1,size(mac,2)-1);
mac_CT = (mac_no_water-mac_only_water)./mac_only_water;

figure(1); clf;
loglog(E*1e3, mac); hold on;
legend(Z);
ylabel('Attenuation in cm^2/g');
xlabel('Photon Energy in KeV');
title({'Photon Attenuation Coefficients for Uranium',...
'see http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z92.html'});

figure(2); clf;
plot(E*1e3, mac_CT); hold on;
legend(Z);
ylabel('CT Number');
xlabel('Photon Energy in KeV');

%%
clear;clc;
Z = {'COMPACT BONE','Blood','Soft Tissue','Water'};
E = [60,100]*1e3;
mac = PhotonAttenuation(Z, log(E), 'mac');

P = PhysProps(Z);

avg_density_vec = cell2mat(P(:,2))';

mac_lin = mac./repmat(avg_density_vec,size(mac,1),1);

% scatter(mac_lin(1,:),mac_lin(2,:))
% legend(Z)

figure;
xlabel('Linear attenuation (\mu) @ E = 60 [KeV]')
ylabel('Linear attenuation (\mu) @ E = 100 [KeV]')

bone_idx = 1;
bone_density_low = 1.8;
bone_density_high = 1.9;
bone_density_diff = bone_density_high - bone_density_low;
bone_pos=[mac_lin(1,bone_idx)*bone_density_low mac_lin(2,bone_idx)*bone_density_low ...
          mac_lin(1,bone_idx)*bone_density_diff mac_lin(2,bone_idx)*bone_density_diff];
rectangle('Position',bone_pos,'FaceColor','r'); hold off;
text(bone_pos(1),bone_pos(2),'Bone','Color','k')

blood_idx = 2;
blood_density_low = 1.056;
blood_density_high = 1.066;
blood_density_diff = blood_density_high - blood_density_low;
blood_pos=[mac_lin(1,blood_idx)*blood_density_low mac_lin(2,blood_idx)*blood_density_low ...
          mac_lin(1,blood_idx)*blood_density_diff mac_lin(2,blood_idx)*blood_density_diff];
rectangle('Position',blood_pos,'FaceColor','b'); hold off;
text(blood_pos(1),blood_pos(2),'Blood','Color','k')

tissue_idx = 3;
tissue_density_low = 0.999;
tissue_density_high = 1.001;
tissue_density_diff = tissue_density_high - tissue_density_low;
tissue_pos=[mac_lin(1,tissue_idx)*tissue_density_low mac_lin(2,tissue_idx)*tissue_density_low ...
          mac_lin(1,tissue_idx)*tissue_density_diff mac_lin(2,tissue_idx)*tissue_density_diff];
rectangle('Position',tissue_pos,'FaceColor','c'); hold off;
text(tissue_pos(1),tissue_pos(2),'Tissue','Color','k')

water_idx = 4;
water_density_low = 0.999;
water_density_high = 1.001;
water_density_diff = water_density_high - water_density_low;
water_pos=[mac_lin(1,water_idx)*water_density_low mac_lin(2,water_idx)*water_density_low ...
          mac_lin(1,water_idx)*water_density_diff mac_lin(2,water_idx)*water_density_diff];
rectangle('Position',water_pos,'FaceColor','g'); hold off;
text(water_pos(1),water_pos(2),'Water','Color','k')

%%
Z = {'COMPACT BONE','Blood','Soft Tissue'};%,'Water'};
E = logspace(log10(0.020), log10(0.150), 500);  % define energy grid
mac = PhotonAttenuation(Z, E, 'mac');

P = PhysProps(Z);

mac_no_water   = mac(:,1:(end-1));
mac_only_water = repmat(mac(:,end),1,size(mac,2)-1);
mac_CT = (mac_no_water-mac_only_water)./mac_only_water;

figure;
loglog(E*1e3, mac,'LineWidth', 5); hold on;
% plot(E*1e3, mac); hold on;
% legend(Z);
set(gca,'FontSize',20);
% ylabel(sprintf('Attenuation\nCoefficients\nin cm^2/g'));
% set(get(gca,'ylabel'),'rotation',0)
% xlabel('Photon Energy in KeV');
% title('Photon Attenuation Coefficients');



