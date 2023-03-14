clear;clc;

energy_low  = 50;
energy_high = 80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Reconstructions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slice_to_plot = 5;
single_init = load(['Single_',num2str(energy_low),'init_vol.mat'],'reconstructed');
reconstruct_low_freq = single_init.reconstructed(:,:,slice_to_plot);
single_init = load(['Single_',num2str(energy_high),'init_vol.mat'],'reconstructed');
reconstruct_high_freq = single_init.reconstructed(:,:,slice_to_plot);
clear single_init;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Reconstructions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,2,1); 
imagesc(reconstruct_low_freq); title('Low Energy - 50 [KeV]');
colormap('bone'); pbaspect([1,1,1]);
subplot(2,2,2); 
imagesc(reconstruct_high_freq); title('High Energy - 80 [KeV]');
colormap('bone'); pbaspect([1,1,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Attenuation Matrix:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gParams;

mu_water_e_low  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_low);
mu_water_e_high = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_high);
mu_bone_e_low   = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_low);
mu_bone_e_high  = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_high);

mu_matrix = [mu_water_e_low,mu_bone_e_low;mu_water_e_high,mu_bone_e_high];

inv_mu_matrix = mu_matrix^-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_mat = inv_mu_matrix*[reconstruct_low_freq(:), reconstruct_high_freq(:)].';
a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3);
imagesc(norm_a1,[0,1]); title('Water Coefficient');
colormap('bone'); pbaspect([1,1,1]);
subplot(2,2,4);
imagesc(norm_a2,[0,1]); title('Bone Coefficient');
colormap('bone'); pbaspect([1,1,1]);