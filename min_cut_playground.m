clear; clc;

load('210920/reconstruction_atten_anti_grid_Enabled.mat');
xcat_gt = load('210920/xcat_reduced.mat');
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

Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater'};

% [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Full_Materials_names_clean);
[linear_atten_mat, ~] = PerBinMaterialsAttenuations(Materials_names);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Attenuation vector:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LoadGParams;

%% Initialization:
clc;
% bin_weights = [1, 1, 1, 1, 1, 1];
% bin_weights = [0, 0.3, 0.6, 0.8, 0.9, 1];
bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];

active_bins = [1,2,3,4,5,6];

active_linear_atten_2d_mat = linear_atten_mat(:,active_bins);
active_bin_weights = bin_weights(active_bins);

reconstructed_bins_vec = zeros(numel(reconstruct),length(active_bins));
for ii=1:length(active_bins)
    reconstructed_bins_vec(:,ii) = reshape(eval(sprintf('reconstruct_bin_%d',active_bins(ii)-1)),[],1);
end

materials_recovered_new = reshape(WeightedNearestNeighbour(active_linear_atten_2d_mat,active_bin_weights,reconstructed_bins_vec),size(reconstruct));

%% Min cut starts here:

single_image_reconstructed = reshape(reconstructed_bins_vec(1:60*60,:),60,60,6);
single_image_material = materials_recovered_new(:,:,1); 
single_image_material_padded = ones(size(single_image_material)+[2,2]);
single_image_material_padded(2:61,2:61) = single_image_material;
single_image_material_neighbours(:,1) = reshape(single_image_material_padded(1:60,2:61),[],1);
single_image_material_neighbours(:,2) = reshape(single_image_material_padded(2:61,1:60),[],1);
single_image_material_neighbours(:,3) = reshape(single_image_material_padded(3:62,2:61),[],1);
single_image_material_neighbours(:,4) = reshape(single_image_material_padded(2:61,3:62),[],1);

query_points_kronned = kron(reshape(single_image_reconstructed,[],6),ones(size(active_linear_atten_2d_mat,1),1));
data_points_repmated = repmat(active_linear_atten_2d_mat,size(reshape(single_image_reconstructed,[],6),1),1);

bin_weights_repmated = repmat(active_bin_weights,size(data_points_repmated,1),1);

diff_per_query = reshape(sqrt(sum(bin_weights_repmated.*((query_points_kronned - data_points_repmated).^2),2)),size(active_linear_atten_2d_mat,1),[]);

factor_per_different_label = 0.03;
label_diff_mat = 4*factor_per_different_label*ones([size(diff_per_query)]);
for ii=1:(60*60)
    for jj=1:8
        label_diff_mat(jj,ii) = label_diff_mat(jj,ii) - factor_per_different_label*(single_image_material_neighbours(ii,1)==jj);
        label_diff_mat(jj,ii) = label_diff_mat(jj,ii) - factor_per_different_label*(single_image_material_neighbours(ii,2)==jj);
        label_diff_mat(jj,ii) = label_diff_mat(jj,ii) - factor_per_different_label*(single_image_material_neighbours(ii,3)==jj);
        label_diff_mat(jj,ii) = label_diff_mat(jj,ii) - factor_per_different_label*(single_image_material_neighbours(ii,4)==jj);
    end
end

[~,nn_old] = min(diff_per_query,[],1);
[~,nn_new] = min(diff_per_query+label_diff_mat,[],1);

figure;
subplot(1,2,1); imagesc(reshape(nn_old,60,60));
subplot(1,2,2); imagesc(reshape(nn_new,60,60));


