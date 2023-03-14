clear; clc;

load('220201\xcat_reduced.mat');

slice = 5;
% xcat_density_slice = xcat_density(:,:,slice);
xcat_id_slice = xcat_id(:,:,slice);

% energy_vec = csvread('spectrum.txt');
% max_energy = max(energy_vec);

% materials_id = unique(xcat_id);
% LoadGParams;
% for ii=1:length(materials_id)
%     Materials_names{ii} = gParamStruct.dict_mat_to_dens_names{find(gParamStruct.dict_mat_to_dens_xcat_id==materials_id(ii))};
%     Materials_names{ii}(Materials_names{ii}==' ') = [];
% end


% %%
load('selected_std.mat');
load('220201\reconstruction_atten_anti_grid_Enabled.mat');
attenuation_slices_noised(:,:,1) = reconstruct_bin_0(:,:,slice);
attenuation_slices_noised(:,:,2) = reconstruct_bin_1(:,:,slice);
attenuation_slices_noised(:,:,3) = reconstruct_bin_2(:,:,slice);
attenuation_slices_noised(:,:,4) = reconstruct_bin_3(:,:,slice);
attenuation_slices_noised(:,:,5) = reconstruct_bin_4(:,:,slice);
attenuation_slices_noised(:,:,6) = reconstruct_bin_5(:,:,slice);
attenuation_full_spectrum_noised = reconstruct(:,:,slice);

% simulation_tissue_std_per_bin = std(reshape(attenuation_slices_noised(repmat(xcat_id_slice,1,1,6)==6),[],6));
% fprintf('Simulation Tissue STD Vec = ');fprintf('%.4f, ',simulation_tissue_std_per_bin);fprintf('\b\b.\n');

%% Recon Method Comparing

LoadGParams;
% bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];
bin_weights = 1./(simulation_tissue_std_per_bin.^1)';
% active_bins = [1,2,3,4,5,6];

% Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater'};
% Full_Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater','Blood','Redmarrow'};
Materials_names = {'Air','Water','Muscle','dryribwater','Redmarrow','Iodineblood','Adiposetissue'};
index_to_id_dict = [23,    1,      2           38             31         18              6];

% [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Full_Materials_names_clean);
[linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names);

energy_vec = csvread('spectrum.txt');
[per_energy_linear_atten_mat, ~] = PerEnergyMaterialsAttenuations(Materials_names, 1:length(energy_vec));
full_spectrum_linear_atten = sum(per_energy_linear_atten_mat(:,:).*repmat(energy_vec',7,1),2);

reconstructed_bins_vec = reshape(attenuation_slices_noised,[],6);

factor_per_different_label = 0.5;
materials_recovered_old = reshape(SmartWeightedNearestNeighbour(full_spectrum_linear_atten,1,attenuation_full_spectrum_noised(:),factor_per_different_label),size(attenuation_slices_noised(:,:,1)));
materials_recovered_mid = reshape(SmartWeightedNearestNeighbour(linear_atten_mat(:,[4,6]),bin_weights([4,6]),reconstructed_bins_vec(:,[4,6]),factor_per_different_label),size(attenuation_slices_noised(:,:,1)));
materials_recovered_new = reshape(SmartWeightedNearestNeighbour(linear_atten_mat,bin_weights,reconstructed_bins_vec,factor_per_different_label),size(attenuation_slices_noised(:,:,1)));

% Compare:
xcat_id_recovered_old = index_to_id_dict(materials_recovered_old);
xcat_id_recovered_mid = index_to_id_dict(materials_recovered_mid);
xcat_id_recovered_new = index_to_id_dict(materials_recovered_new);

figure;
% subplot(2,2,1);imagesc(xcat_id_recovered_old); pbaspect([1,1,1]); title('Full Spectrum Detectors');
% subplot(2,2,2);imagesc(xcat_id_recovered_mid); pbaspect([1,1,1]); title('2 Energy Bins from 6');
% subplot(2,2,3);imagesc(xcat_id_recovered_new); pbaspect([1,1,1]); title('6 Energy Bins from 6');
% subplot(2,2,4);imagesc(xcat_id_slice); pbaspect([1,1,1]); title('Ground Truth');

subplot(1,3,1);
imagesc(xcat_id_recovered_old, [1, 38]);
pbaspect([1,1,1]); title('Full Spectrum Detectors');
subplot(1,3,2);
imagesc(xcat_id_recovered_new, [1, 38]);
pbaspect([1,1,1]); title('Multi Spectral Detectors');
subplot(1,3,3);
imagesc(xcat_id_slice, [1, 38]);
pbaspect([1,1,1]); title('Ground Truth');

%%

figure;
subplot(1,2,1);
imshow(attenuation_full_spectrum_noised,[]);
title('Full Spectrum Attenuation Map');
for ii=1:6
    subplot(2,6,ii+3+3*(ii>3));
    imshow(attenuation_slices_noised(:,:,ii),[]);
    title(sprintf('Bin %d',ii));
end

