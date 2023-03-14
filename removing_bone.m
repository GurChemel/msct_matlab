clear; clc;

load('reconstruction_atten_per_bin1.mat')

%%

original_image = reconstruct_bin_3;
final_image = reconstruct_bin_3;
bone_image = zeros(size(final_image));

figure;
subplot(1,2,1);
imagesc(squeeze(sum(permute(original_image,[1,2,3]),3)));
title('Original Image');
colormap('bone'); pbaspect([1,1,1]);

%%

selected_bins = [4,6];

energy_centers = [31.3090   41.5364   50.3864   58.9392   69.0782   90.8517];

energy_low  = energy_centers(selected_bins(1));
energy_high = energy_centers(selected_bins(2));

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

for slice = 1:10
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loading Reconstructions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    reconstruct_low_freq  = reconstruct_bin_3(:,:,slice);
    reconstruct_high_freq = reconstruct_bin_5(:,:,slice);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting Decompositions:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_mat = inv_mu_matrix*[reconstruct_low_freq(:), reconstruct_high_freq(:)].';
    a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
    a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

    norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
    norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);

    bone_idx = (0.8*norm_a1)<norm_a2;
    bone_idx = medfilt2(bone_idx, [2 2]);

    final_image(:,:,slice) = (1-bone_idx).*original_image(:,:,slice);
    bone_image(:,:,slice) = bone_idx;
    
end

%%
subplot(1,2,2);
imagesc(squeeze(sum(permute(final_image,[1,2,3]),3)));
title('Boneless Image');
colormap('bone'); pbaspect([1,1,1]);

%%
isosurface(bone_image,1)