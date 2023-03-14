clear; clc;

load('asg_reconstruction_atten_low_high.mat');

reconstruct_high_freq = asg_reconstruct_high_freq(:,:,5);
reconstruct_low_freq  = asg_reconstruct_low_freq(:,:,5);

gParams;
%%
load('ElementXS.mat');

energy_thresh = 38; % KeV
energy_vec = csvread('spectrum.txt');
% energy_vec(150)=[];
energy_vec(150)=[];

Na = 6.02214076e23;

% idx       =   1   2   3   4   5   6   7
% Names_vec = {'H','C','N','O','P','K','Ca'};
% wanted_z  = [ 1 , 6 , 7 , 8 , 15, 19, 20];

zvec_h_idx  = find(ZVec==1,1);
zvec_o_idx  = find(ZVec==8,1);
zvec_p_idx  = find(ZVec==15,1);
zvec_ca_idx = find(ZVec==20,1);

A_h = 1.008;
A_o = 15.999;
A_p = 30.974;
A_ca = 40.078;

h_xs_per_e  = Photon(:,zvec_h_idx).*energy_vec*Na/A_h;
o_xs_per_e  = Photon(:,zvec_o_idx).*energy_vec*Na/A_o;
p_xs_per_e  = Photon(:,zvec_p_idx).*energy_vec*Na/A_p;
ca_xs_per_e = Photon(:,zvec_ca_idx).*energy_vec*Na/A_ca;

mu_water_e_low_calculated  = 2*sum(h_xs_per_e(1:energy_thresh))+sum(o_xs_per_e(1:energy_thresh));
mu_water_e_high_calculated = 2*sum(h_xs_per_e((energy_thresh):end))+sum(o_xs_per_e((energy_thresh):end));

mu_bone_e_low_calculated  = (sum(ca_xs_per_e(1:energy_thresh))+0.58*sum(p_xs_per_e(1:energy_thresh)));
mu_bone_e_high_calculated = (sum(ca_xs_per_e((energy_thresh):end))+0.58*sum(p_xs_per_e((energy_thresh):end)));

%%

energy_thresh = 60;
max_energy = length(energy_vec);

energy_low  = sum(energy_vec(1:energy_thresh)'.*(1:energy_thresh))/sum(energy_vec(1:energy_thresh));           % In KeV;
energy_high = sum(energy_vec(energy_thresh:max_energy)'.*(energy_thresh:max_energy))/sum(energy_vec(energy_thresh:max_energy)); % In KeV;

% energy_low  = 44;
% energy_high = 77;

% water_mass_density_in_room_temp = 0.99802;
% water_atten(:,2) = water_mass_density_in_room_temp*water_atten(:,2);
% bone_mass_density_in_room_temp = 1.8;
% bone_atten(:,2)  = bone_mass_density_in_room_temp*bone_atten(:,2); % bone_Attem = mu/ro

mu_water_e_low  = sum(interp1((1e3)*water_atten(:,1),water_atten(:,3),1:energy_thresh).*energy_vec(1:energy_thresh)')/sum(energy_vec(1:energy_thresh));
mu_water_e_high = sum(interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_thresh:max_energy).*energy_vec(energy_thresh:max_energy)')/sum(energy_vec(energy_thresh:max_energy));
mu_bone_e_low   = sum(interp1((1e3)*bone_atten(:,1),bone_atten(:,3),1:energy_thresh).*energy_vec(1:energy_thresh)')/sum(energy_vec(1:energy_thresh));
mu_bone_e_high  = sum(interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_thresh:max_energy).*energy_vec(energy_thresh:max_energy)')/sum(energy_vec(energy_thresh:max_energy));

% mu_water_e_low  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_low);
% mu_water_e_high = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_high);
% mu_bone_e_low   = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_low);
% mu_bone_e_high  = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_high);

mu_matrix = [mu_water_e_low,mu_bone_e_low;mu_water_e_high,mu_bone_e_high];
% mu_matrix = [mu_water_e_low_calculated,mu_bone_e_low_calculated;mu_water_e_high_calculated,mu_bone_e_high_calculated];

inv_mu_matrix = mu_matrix^-1;

%%

a_mat = inv_mu_matrix*[reconstruct_low_freq(:), reconstruct_high_freq(:)].';
a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

% a_1 = inv_mu_matrix(1,1)*reconstruct_low_freq+inv_mu_matrix(1,2)*reconstruct_high_freq;
% a_2 = inv_mu_matrix(2,1)*reconstruct_low_freq+inv_mu_matrix(2,2)*reconstruct_high_freq;

%%
% loglog(ct_to_dens(:,1),ct_to_dens(:,2))

% val=2;
% 
% figure;
% subplot(1,2,1); imshow(conv2(reconstruct_high_freq,(1/(val^2))*ones(val,val)),[0,0.6]); title('High Freq');
% subplot(1,2,2); imshow(conv2(reconstruct_low_freq,(1/(val^2))*ones(val,val)),[0,0.6]); title('Low Freq');

%%
norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);

% % vector_of_xy_values = (1:50) - 25;
% % [Yg, Xg] = ndgrid(vector_of_xy_values, vector_of_xy_values);
% % idx = find(Xg.^2 + Yg.^2 > 180);
% % 
% % norm_a1(idx)=0;
% % norm_a2(idx)=0;
% 
% figure;
% subplot(2,2,1); imshow(reconstruct_high_freq,[0,0.6]); title('High Freq');
% subplot(2,2,3); imshow(reconstruct_low_freq ,[0,0.6]); title('Low Freq');
% subplot(2,2,2); imshow(norm_a1,[0,1]); title('Water');
% subplot(2,2,4); imshow(norm_a2,[0,1]); title('Bone');

%%
% figure;
% subplot(2,3,1); imshow(reconstruct_high_freq,[0,0.6]); title('High Freq');
% subplot(2,3,4); imshow(reconstruct_low_freq ,[0,0.6]); title('Low Freq');
% % subplot(2,3,3); imshow(a_1,[0,max(a_1,[],'all')]); title('Water');
% % subplot(2,3,6); imshow(a_2,[0,max(a_2,[],'all')]); title('Bone');
% subplot(2,3,3); imshow(norm_a1,[0,1]); title('Water');
% subplot(2,3,6); imshow(norm_a2,[0,1]); title('Bone');
% subplot(2,3,2); imshow(a_1*mu_water_e_high+a_2*mu_bone_e_high,[0,0.6]); title('Rec High Freq');
% subplot(2,3,5); imshow(a_1*mu_water_e_low+a_2*mu_bone_e_low,[0,0.6]); title('Rec Low Freq');
% % subplot(2,3,2); imshow(reconstruct_high_freq - a_1*mu_water_e_high+a_2*mu_bone_e_high,[0,0.6]); title('Error High Freq');
% % subplot(2,3,5); imshow(reconstruct_low_freq - a_1*mu_water_e_low+a_2*mu_bone_e_low,[0,0.6]); title('Error Low Freq');

%%
figure;
subplot(2,2,1); imshow(reconstruct_high_freq,[0,0.6]); title('High Energy Reconstruction');
subplot(2,2,3); imshow(reconstruct_low_freq ,[0,0.6]); title('Low Energy Reconstruction');
subplot(2,2,2); imshow(norm_a1,[0,1]); title('Water Coefficient');
subplot(2,2,4); imshow(norm_a2,[0,1]); title('Bone Coefficient');
%%
% figure;
% subplot(2,3,1);
% imagesc(reconstruct_high_freq,[0,0.6]);
% title('High Energy Reconstruction');
% pbaspect([1,1,1]);
% set(gca,'xtick',[]); set(gca,'ytick',[]);
% 
% subplot(2,3,4);
% imagesc(reconstruct_low_freq ,[0,0.6]);
% title('Low Energy Reconstruction');
% pbaspect([1,1,1]);
% set(gca,'xtick',[]); set(gca,'ytick',[]);
% 
% subplot(2,3,2);
% imagesc(norm_a1);%,[0,1]);
% title('Water Coefficient');
% pbaspect([1,1,1]);
% set(gca,'xtick',[]); set(gca,'ytick',[]);
% 
% subplot(2,3,5);
% imagesc(norm_a2);%,[0,1]);
% title('Bone Coefficient');
% pbaspect([1,1,1]);
% set(gca,'xtick',[]); set(gca,'ytick',[]);
% 
% subplot(2,3,3);
% xcat_water = zeros(size(xcat_density));
% xcat_water(((xcat_id==1) + (xcat_id==2))>0) = xcat_density(((xcat_id==1) + (xcat_id==2))>0);
% imagesc(xcat_water);%,[0,1]);
% title('GT - Water Coefficient');
% pbaspect([1,1,1]);
% set(gca,'xtick',[]); set(gca,'ytick',[]);
% 
% subplot(2,3,6);
% xcat_bone = zeros(size(xcat_density));
% xcat_bone(xcat_id==3) = xcat_density(xcat_id==3);
% imagesc(xcat_bone);%,[0,1]);
% title('GT - Bone Coefficient');
% pbaspect([1,1,1]);
% set(gca,'xtick',[]); set(gca,'ytick',[]);
