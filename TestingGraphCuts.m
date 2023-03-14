clear; clc;
addpath('GCMex/')

% load('220201\xcat_reduced.mat');
% 
% slice = 5;
% xcat_density_slice = xcat_density(:,:,slice);
% xcat_id_slice = xcat_id(:,:,slice);
% 
% energy_vec = csvread('spectrum.txt');
% max_energy = max(energy_vec);
% 
% materials_id = unique(xcat_id);
% LoadGParams;
% for ii=1:length(materials_id)
%     Materials_names{ii} = gParamStruct.dict_mat_to_dens_names{find(gParamStruct.dict_mat_to_dens_xcat_id==materials_id(ii))};
%     Materials_names{ii}(Materials_names{ii}==' ') = [];
% end
% 
% % [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Full_Materials_names_clean);
% [per_bin_linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names);
% [per_energy_linear_atten_mat, ~] = PerEnergyMaterialsAttenuations(Materials_names, 1:length(energy_vec));
% 
% attenuation_slices = zeros([size(xcat_density_slice),6]);
% attenuation_full_spectrum = zeros([size(xcat_density_slice)]);
% for row = 1:size(xcat_density_slice,1)
%     for col = 1:size(xcat_density_slice,2)
% %         attenuation_slices(row,col,:) = xcat_density_slice(row,col)*per_bin_linear_atten_mat(materials_id==xcat_id_slice(row,col),:);
% %         attenuation_full_spectrum(row,col) = xcat_density_slice(row,col)*sum(per_energy_linear_atten_mat(materials_id==xcat_id_slice(row,col),:).*energy_vec');
%         attenuation_slices(row,col,:) = per_bin_linear_atten_mat(materials_id==xcat_id_slice(row,col),:);
%         attenuation_full_spectrum(row,col) = sum(per_energy_linear_atten_mat(materials_id==xcat_id_slice(row,col),:).*energy_vec');
%     end
% end

%%
% clc;
% load('selected_std.mat');
% 
% % load('220201\reconstruction_atten_anti_grid_Enabled.mat');
% % attenuation_slices_noised(:,:,1) = reconstruct_bin_0(:,:,slice);
% % attenuation_slices_noised(:,:,2) = reconstruct_bin_1(:,:,slice);
% % attenuation_slices_noised(:,:,3) = reconstruct_bin_2(:,:,slice);
% % attenuation_slices_noised(:,:,4) = reconstruct_bin_3(:,:,slice);
% % attenuation_slices_noised(:,:,5) = reconstruct_bin_4(:,:,slice);
% % attenuation_slices_noised(:,:,6) = reconstruct_bin_5(:,:,slice);
% 
% noise_level = 0.25;
% noise_factor = permute([1,1.1,1.35,1.5,1.8,2],[1,3,2]);
% attenuation_slices_noised = attenuation_slices + noise_level*repmat(noise_factor.*max(attenuation_slices,[],[1,2]),60,60,1).*(rand(size(attenuation_slices))-0.5);
% 
% fake_std_vec = std(reshape(attenuation_slices_noised(repmat(xcat_id_slice,1,1,6)==6),[],6));
% fprintf('Simulation Tissue STD Vec = ');fprintf('%.4f, ',simulation_tissue_std_per_bin);fprintf('\b\b.\n');
% fprintf('Fake Image Tissue STD Vec = ');fprintf('%.4f, ',fake_std_vec);fprintf('\b\b.\n');

%% Recon Method Comparing
% 
% addpath('GCMex/')
% 
% LoadGParams;
% bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];
% 
% % active_bins = [1,2,3,4,5,6];
% 
% % Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater'};
% % Full_Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater','Blood','Redmarrow'};
% Materials_names = {'Air','Water','Muscle','dryribwater','Redmarrow','Iodineblood','Adiposetissue'};
% index_to_id_dict = [23,    1,      2           38             31         18              6];
% 
% % [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Full_Materials_names_clean);
% [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names);
% 
% energy_vec = csvread('spectrum.txt');
% [per_energy_linear_atten_mat, ~] = PerEnergyMaterialsAttenuations(Materials_names, 1:length(energy_vec));
% full_spectrum_linear_atten = sum(per_energy_linear_atten_mat(:,:).*repmat(energy_vec',7,1),2);
% 
% % active_linear_atten_2d_mat = linear_atten_mat(:,active_bins);
% % active_bin_weights = bin_weights(active_bins);
% 
% % noise_level_vec = 0.2:0.1:0.4;
% % noise_level_vec = 0:0.01:0.2;
% noise_level_vec = 0:0.4:0.8;
% 
% num_runs_per_noise = 10;
% 
% error_old_vec = zeros(num_runs_per_noise,length(noise_level_vec));
% error_mid_vec = zeros(num_runs_per_noise,length(noise_level_vec));
% error_new_vec = zeros(num_runs_per_noise,length(noise_level_vec));
% error_gct_vec = zeros(num_runs_per_noise,length(noise_level_vec));
% 
% psnr_vec = zeros(num_runs_per_noise,length(noise_level_vec));
% 
% old_recon_database = zeros(60,60,length(noise_level_vec));
% mid_recon_database = zeros(60,60,length(noise_level_vec));
% new_recon_database = zeros(60,60,length(noise_level_vec));
% gct_recon_database = zeros(60,60,length(noise_level_vec));
% 
% for run_num = 1:num_runs_per_noise
%     for jj=1:length(noise_level_vec)
% 
%         noise_level = noise_level_vec(jj);
% 
%         attenuation_full_spectrum_noised = attenuation_full_spectrum + noise_level*max(attenuation_full_spectrum,[],[1,2])*(rand(size(attenuation_full_spectrum))-0.5);
%         attenuation_slices_noised = attenuation_slices + noise_level*repmat(max(attenuation_slices,[],[1,2]),60,60,1).*(rand(size(attenuation_slices))-0.5);
% 
%         reconstructed_bins_vec = reshape(attenuation_slices_noised,[],6);
% 
%         % materials_recovered_new = reshape(knnsearch(active_linear_atten_2d_mat,reconstructed_bins_vec,'Distance','cityblock'),size(attenuation_slices_noised(:,:,1)));
%         % materials_recovered_new = reshape(WeightedNearestNeighbour(active_linear_atten_2d_mat,active_bin_weights,reconstructed_bins_vec),size(attenuation_slices_noised(:,:,1)));
% 
%         factor_per_different_label = 0.5;
%         smothness_factor = 0.03;
% %         materials_recovered_old = reshape(SmartWeightedNearestNeighbour(full_spectrum_linear_atten,1,attenuation_full_spectrum_noised(:),factor_per_different_label),size(attenuation_slices_noised(:,:,1)));
%         materials_recovered_old = GraphCutAccumWrapper(attenuation_full_spectrum_noised, smothness_factor);
%         materials_recovered_mid = reshape(SmartWeightedNearestNeighbour(linear_atten_mat(:,[4,6]),bin_weights([4,6]),reconstructed_bins_vec(:,[4,6]),factor_per_different_label),size(attenuation_slices_noised(:,:,1)));
%         materials_recovered_new = reshape(SmartWeightedNearestNeighbour(linear_atten_mat,bin_weights,reconstructed_bins_vec,factor_per_different_label),size(attenuation_slices_noised(:,:,1)));
%         materials_recovered_gct = GraphCutPerBinWrapper(attenuation_slices_noised, smothness_factor);
%         
%         % Compare:
%         xcat_id_recovered_old = index_to_id_dict(materials_recovered_old);
%         xcat_id_recovered_mid = index_to_id_dict(materials_recovered_mid);
%         xcat_id_recovered_new = index_to_id_dict(materials_recovered_new);
%         xcat_id_recovered_gct = index_to_id_dict(materials_recovered_gct);
% 
%         if run_num==1
%             old_recon_database(:,:,jj) = xcat_id_recovered_old;
%             mid_recon_database(:,:,jj) = xcat_id_recovered_mid;
%             new_recon_database(:,:,jj) = xcat_id_recovered_new;
%             gct_recon_database(:,:,jj) = xcat_id_recovered_gct;
%         end
%         
%         psnr_vec(run_num,jj) = psnr(attenuation_full_spectrum_noised,attenuation_full_spectrum);
%         
%         error_old_vec(run_num,jj) = sum(1-(xcat_id_recovered_old==xcat_id_slice),'all')/numel(xcat_id_slice);
%         error_mid_vec(run_num,jj) = sum(1-(xcat_id_recovered_mid==xcat_id_slice),'all')/numel(xcat_id_slice);
%         error_new_vec(run_num,jj) = sum(1-(xcat_id_recovered_new==xcat_id_slice),'all')/numel(xcat_id_slice);
%         error_gct_vec(run_num,jj) = sum(1-(xcat_id_recovered_gct==xcat_id_slice),'all')/numel(xcat_id_slice);
%         
%         if run_num == 1
%             if jj == 1
%                 figure;
%             end
%             subplot(length(noise_level_vec),4,4*(jj-1)+1); imagesc(xcat_id_recovered_old); pbaspect([1,1,1]); title('Full Spectrum Detectors');
%             subplot(length(noise_level_vec),4,4*(jj-1)+2); imagesc(xcat_id_recovered_new); pbaspect([1,1,1]); title('Multi Spectral Detectors');
%             subplot(length(noise_level_vec),4,4*(jj-1)+3); imagesc(xcat_id_recovered_gct); pbaspect([1,1,1]); title('Multi Spectral Detectors - GraphCut');
%             subplot(length(noise_level_vec),4,4*(jj-1)+4); imagesc(xcat_id_slice); pbaspect([1,1,1]); title('Ground Truth');
%         end
% 
%     end
% end

%% Recon Method Comparing

% load('220426\reconstruction_atten_anti_grid_Disabled.mat')

load('220621\reconstruction_full_spectrum.mat')
load('220621\reconstruction_atten_per_bin0.mat')
reconstruct = reconstruct_with_asg;
reconstruct_per_bin(:,:,:,1) = reconstruct_bin_0;
reconstruct_per_bin(:,:,:,2) = reconstruct_bin_1;
reconstruct_per_bin(:,:,:,3) = reconstruct_bin_2;
reconstruct_per_bin(:,:,:,4) = reconstruct_bin_3;
reconstruct_per_bin(:,:,:,5) = reconstruct_bin_4;
reconstruct_per_bin(:,:,:,6) = reconstruct_bin_5;

% load('220510\reconstruction_150s.mat')
% reconstruct = X_reshaped*2.8;
% load('220517\reconstruction_per_bin_150s.mat')
% reconstruct_per_bin = X_reshaped*2.8;

LoadGParams;
bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];

% Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater'};
% Full_Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater','Blood','Redmarrow'};
% Materials_names = {'Air','Water','Muscle','dryribwater','Redmarrow','Iodineblood','Adiposetissue'};
index_to_id_dict = [23,    1,      2           38             31         18              6];
% [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names);
% 
% energy_vec = csvread('spectrum.txt');
% [per_energy_linear_atten_mat, ~] = PerEnergyMaterialsAttenuations(Materials_names, 1:length(energy_vec));
% full_spectrum_linear_atten = sum(per_energy_linear_atten_mat(:,:).*repmat(energy_vec',7,1),2);

smothness_factor_vec = 0.00:0.01:0.04;

for jj=1:length(smothness_factor_vec)

    factor_per_different_label = 0.5;
    smothness_factor = smothness_factor_vec(jj);
    materials_recovered_old = GraphCutAccumWrapper(squeeze(reconstruct(:,:,8)), smothness_factor);
    materials_recovered_gct = GraphCutPerBinWrapper(squeeze(reconstruct_per_bin(:,:,8,:)), smothness_factor);

    % Compare:
    xcat_id_recovered_old = index_to_id_dict(materials_recovered_old);
    xcat_id_recovered_gct = index_to_id_dict(materials_recovered_gct);

    if jj == 1
        figure;
    end
    subplot(2,length(smothness_factor_vec),jj); imagesc(xcat_id_recovered_old); pbaspect([1,1,1]); title('Full Spectrum Detectors');
    subplot(2,length(smothness_factor_vec),jj+length(smothness_factor_vec)); imagesc(xcat_id_recovered_gct); pbaspect([1,1,1]); title('Multi Spectral Detectors - GraphCut');
end
