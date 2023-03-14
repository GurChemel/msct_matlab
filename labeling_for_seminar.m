clear; clc;


load('221117/reconstruction_atten_per_bin0.mat');
slice = 5;

image_values = zeros([size(reconstruct_bin_0,1),size(reconstruct_bin_0,2),6]);
figure;
for bin=1:6
    image_values(:,:,bin) = eval(sprintf('reconstruct_bin_%d(:,:,%d)',bin-1,slice));
end
image_values = image_values/max(image_values,[],'all');
image_values = (image_values.^0.5);
for bin=1:6
    subplot(2,3,bin);
    imshow(image_values(:,:,bin));
end

attenuation_slices_noised(:,:,1) = reconstruct_bin_0(:,:,slice);
attenuation_slices_noised(:,:,2) = reconstruct_bin_1(:,:,slice);
attenuation_slices_noised(:,:,3) = reconstruct_bin_2(:,:,slice);
attenuation_slices_noised(:,:,4) = reconstruct_bin_3(:,:,slice);
attenuation_slices_noised(:,:,5) = reconstruct_bin_4(:,:,slice);
attenuation_slices_noised(:,:,6) = reconstruct_bin_5(:,:,slice);

%% Recon Method Comparing


addpath('GCMex/')

LoadGParams;

% bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];
% 
% % active_bins = [1,2,3,4,5,6];
% 
% Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater'};
% Full_Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater','Blood','Redmarrow'};
Materials_names = {'Air','Water','Muscle','dryribwater','Redmarrow','Iodineblood','Adiposetissue'};
index_to_id_dict = [23,    1,      2           38             31         18              6];
% 
% % [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Full_Materials_names_clean);
% [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names);
% 
% energy_vec = csvread('spectrum.txt');
% [per_energy_linear_atten_mat, ~] = PerEnergyMaterialsAttenuations(Materials_names, 1:length(energy_vec));
% full_spectrum_linear_atten = sum(per_energy_linear_atten_mat(:,:).*repmat(energy_vec',7,1),2);

smothness_factor = 0.01;
materials_recovered_gct = GraphCutPerBinWrapper(attenuation_slices_noised, smothness_factor);

xcat_id_recovered_gct = index_to_id_dict(materials_recovered_gct);

imagesc(xcat_id_recovered_gct)
