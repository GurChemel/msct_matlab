clear;clc;

%%
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
mono_a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
mono_a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

%%

load('reconstruction_atten_per_bin0.mat')

per_bin_image(:,:,1) = reconstruct_bin_0(:,:,slice_to_plot);
per_bin_image(:,:,2) = reconstruct_bin_1(:,:,slice_to_plot);
per_bin_image(:,:,3) = reconstruct_bin_2(:,:,slice_to_plot);
per_bin_image(:,:,4) = reconstruct_bin_3(:,:,slice_to_plot);
per_bin_image(:,:,5) = reconstruct_bin_4(:,:,slice_to_plot);
per_bin_image(:,:,6) = reconstruct_bin_5(:,:,slice_to_plot);

selected_bins = [4,6];

energy_centers = [31.3090   41.5364   50.3864   58.9392   69.0782   90.8517];

energy_low  = energy_centers(selected_bins(1));
energy_high = energy_centers(selected_bins(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Reconstructions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reconstruct_low_freq  = per_bin_image(:,:,selected_bins(1));
reconstruct_high_freq = per_bin_image(:,:,selected_bins(2));

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
multi_a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
multi_a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

%%
load('xcat_reduced.mat')

slice_to_plot = 5;
xcat_id_slice = xcat_id(:,:,slice_to_plot);
bone_gt_id = (xcat_id_slice==38);
water_gt_id = ((xcat_id_slice~=23).*(xcat_id_slice~=38))==1;

% gt_bone_density_vec = 

% subplot(1,2,1); imshow(bone_gt_id.*xcat_density(:,:,5))
% subplot(1,2,2); imshow(water_gt_id.*xcat_density(:,:,5))

% write_im_for_pres(bone_gt_id,'bone_id_gt')
% write_im_for_pres(water_gt_id,'water_gt_id')
