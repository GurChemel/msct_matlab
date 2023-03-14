clear; clc;

load('220113/reconstruction_atten_anti_grid_Enabled.mat');
xcat_gt = load('220113/xcat_reduced.mat');
materials_db_cell = load('materials.mat');
materials_id = load('xcat_materials.mat');
Full_Materials_names = fields(materials_id)';
Full_Materials_names_clean = Full_Materials_names;
materials_db_mat = zeros(length(Full_Materials_names),length(materials_db_cell.Air));
for ii=1:length(Full_Materials_names)
    Full_Materials_names_clean{ii}(Full_Materials_names{ii}=='_') = '';
    materials_db_mat(ii,:) = materials_db_cell.(Full_Materials_names_clean{ii});
end
materials_db_mat(:,end) = [];

% for ii=1:6
%     write_im_for_pres(eval(sprintf('reconstruct_bin_%d(:,:,1)',ii-1)),sprintf('smart_weight_knn/reconstruct_bin_%d',ii-1),4)
% end

% Reduce Materials:
% active_materials = [1,8,14,19];
% materials_to_delete = 1:27;
% materials_to_delete(active_materials)=[];
% materials_db_mat(:,materials_to_delete)=0;
% materials_db_mat = materials_db_mat./sum(materials_db_mat,2);

% Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater'};
% Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater','Blood','Redmarrow'};
Materials_names = {'Air','Water','Muscle','dryribwater','Redmarrow'};

% [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Full_Materials_names_clean);
[linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names);

% linear_atten_mat = linear_atten_mat.*repmat(density_vec',1,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Attenuation vector:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LoadGParams;
full_center_energy = 50;
full_mu_water_e     = interp1((1e3)*gParamStruct.water_atten(:,1),0.99857*gParamStruct.water_atten(:,2),full_center_energy);
full_recovery_ct = 1e3*((reconstruct-full_mu_water_e)/(full_mu_water_e));
full_recovery_density = interp1(gParamStruct.ct_to_dens(:,1),gParamStruct.ct_to_dens(:,2),full_recovery_ct);

%% Old Method:
materials_recovered_old = zeros(size(full_recovery_density));
for ii=1:size(gParamStruct.dict_mat_to_dens_values,1)
    materials_recovered_old = materials_recovered_old+ii*(full_recovery_density>=gParamStruct.dict_mat_to_dens_values(ii,1)).*(full_recovery_density<gParamStruct.dict_mat_to_dens_values(ii,2));
end

% Only air-water-bone:
materials_recovered_old_awb = 23*ones(size(full_recovery_density));
materials_recovered_old_awb(full_recovery_density>=gParamStruct.dict_mat_to_dens_values(1,2)) = 1;
materials_recovered_old_awb(full_recovery_density>=gParamStruct.dict_mat_to_dens_values(8,1)) = 26;

%% New Method:
clc;
% bin_weights = [1, 1, 1, 1, 1, 1];
% bin_weights = [0, 0, 0, 0.8, 0.9, 1];

bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];

active_bins = [1,2,3,4,5,6];

active_linear_atten_2d_mat = linear_atten_mat(:,active_bins);
active_bin_weights = bin_weights(active_bins);

reconstructed_bins_vec = zeros(numel(reconstruct),length(active_bins));
for ii=1:length(active_bins)
    reconstructed_bins_vec(:,ii) = reshape(eval(sprintf('reconstruct_bin_%d',active_bins(ii)-1)),[],1);
end

% materials_recovered_new = reshape(knnsearch(active_linear_atten_2d_mat,reconstructed_bins_vec,'Distance','cityblock'),size(reconstruct));
% materials_recovered_new = reshape(WeightedNearestNeighbour(active_linear_atten_2d_mat,active_bin_weights,reconstructed_bins_vec),size(reconstruct));

factor_per_different_label = 0.5;
materials_recovered_new = reshape(SmartWeightedNearestNeighbour(active_linear_atten_2d_mat,active_bin_weights,reconstructed_bins_vec,factor_per_different_label),size(reconstruct));

% Compare:
xcat_id_recovered_old = gParamStruct.dict_mat_to_dens_xcat_id(materials_recovered_old);
xcat_id_recovered_new = gParamStruct.dict_mat_to_dens_xcat_id(materials_recovered_new);

error_old = sum(1-(xcat_id_recovered_old==xcat_gt.xcat_id),'all')/numel(xcat_gt.xcat_id);
error_new = sum(1-(xcat_id_recovered_new==xcat_gt.xcat_id),'all')/numel(xcat_gt.xcat_id);

fprintf('Old ID Errors: %d%%\n',round(error_old*100));
fprintf('New ID Errors: %d%%\n',round(error_new*100));

fractional_recovered_old = materials_db_mat(xcat_id_recovered_old(:),:);
fractional_densities_recovered_old = fractional_recovered_old.*repmat(full_recovery_density(:),1,size(fractional_recovered_old,2));

fractional_recovered_new = materials_db_mat(xcat_id_recovered_new(:),:);
fractional_densities_recovered_new = fractional_recovered_new.*repmat(full_recovery_density(:),1,size(fractional_recovered_new,2));

fractional_densities_xcat_gt = reshape(permute(xcat_gt.xcat_fractional_mass,[2,3,4,1]),[],27).*repmat(xcat_gt.xcat_density(:),1,size(fractional_recovered_new,2));

active_materials = 1:27; %[1,8,14,19];

density_error_old = sum((fractional_densities_recovered_old(:,active_materials)-fractional_densities_xcat_gt(:,active_materials)).^2,'all');
density_error_new = sum((fractional_densities_recovered_new(:,active_materials)-fractional_densities_xcat_gt(:,active_materials)).^2,'all');

fprintf('Old Density Errors: %d\n',density_error_old);
fprintf('New Density Errors: %d\n',density_error_new);
fprintf('Density Errors new/old: %d%%\n',round(100*density_error_new/density_error_old));

%%

old_image  = reshape(fractional_densities_recovered_old(1:(60*60),:),60,60,27);
new_image  = reshape(fractional_densities_recovered_new(1:(60*60),:),60,60,27);
xcat_image = reshape(fractional_densities_xcat_gt(1:(60*60),:),60,60,27);

% active_materials = [1,8,14,19];

old_image_normed = zeros(size(old_image));
old_image_normed(:,:,active_materials) = old_image(:,:,active_materials)./sum(old_image(:,:,active_materials),3);

new_image_normed = zeros(size(new_image));
new_image_normed(:,:,active_materials) = new_image(:,:,active_materials)./sum(new_image(:,:,active_materials),3);

xcat_image_normed = zeros(size(xcat_image));
xcat_image_normed(:,:,active_materials) = xcat_image(:,:,active_materials)./sum(xcat_image(:,:,active_materials),3);

max_value_vec = double(max([max(fractional_densities_recovered_old(1:(60*60),:)); ...
                            max(fractional_densities_recovered_new(1:(60*60),:)); ...
                            max(fractional_densities_xcat_gt(1:(60*60),:))]));

max_value_vec(max_value_vec==0) = 1;

num_active_materials = length(active_materials);

elements_to_plot = find(squeeze(sum(sum(xcat_image,1),2)));

figure;
for ii=1:length(elements_to_plot)
    mat_idx = elements_to_plot(ii);
    subplot(3,length(elements_to_plot),ii)
    imagesc(old_image(:,:,mat_idx),[0, max_value_vec(mat_idx)]);
    if ii==1
        title(sprintf('OLD %d : %.2d',mat_idx,max_value_vec(mat_idx)));
    else
        title(sprintf('%d : %.2d',mat_idx,max_value_vec(mat_idx)));
    end
    subplot(3,length(elements_to_plot),length(elements_to_plot)+ii)
    imagesc(new_image(:,:,mat_idx),[0, max_value_vec(mat_idx)]);
    if ii==1
        title('New');
    end
    subplot(3,length(elements_to_plot),2*length(elements_to_plot)+ii)
    imagesc(xcat_image(:,:,mat_idx),[0, max_value_vec(mat_idx)]);
    if ii==1
        title('Xcat');
    end
end


%%
figure;
subplot(2,2,1); imagesc(xcat_id_recovered_old(:,:,1));
subplot(2,2,2); imagesc(xcat_id_recovered_new(:,:,1));
subplot(2,2,3); imagesc(xcat_id_recovered_old(:,:,1));
subplot(2,2,4); imagesc(xcat_id_recovered_new(:,:,1));
%%
figure;
subplot(2,2,1); imagesc(full_recovery_density(:,:,1).*(xcat_id_recovered_new(:,:,1)==26)); colormap('bone'); title('new bone');
subplot(2,2,2); imagesc(full_recovery_density(:,:,1).*((xcat_id_recovered_new(:,:,1)==3) + (xcat_id_recovered_new(:,:,1)==6))); colormap('bone'); title('new water');
subplot(2,2,3); imagesc(full_recovery_density(:,:,1).*(materials_recovered_old_awb(:,:,1)==26)); colormap('bone'); title('old bone');
subplot(2,2,4); imagesc(full_recovery_density(:,:,1).*(materials_recovered_old_awb(:,:,1)==1)); colormap('bone'); title('old water');
%%
% write_im_for_pres(full_recovery_density(:,:,1).*(xcat_id_recovered_new(:,:,1)==26),'smart_weight_knn/new_bone_coeff',4)
% write_im_for_pres(full_recovery_density(:,:,1).*((xcat_id_recovered_new(:,:,1)==3) + (xcat_id_recovered_new(:,:,1)==6)),'smart_weight_knn/new_water_coeff',4)
% 
% write_im_for_pres(full_recovery_density(:,:,1).*(xcat_id_recovered_old(:,:,1)==38),'smart_weight_knn/old_bone_coeff',4)
% write_im_for_pres(full_recovery_density(:,:,1).*((xcat_id_recovered_old(:,:,1)~=38).*(xcat_id_recovered_old(:,:,1)~=1)),'smart_weight_knn/old_water_coeff',4)
