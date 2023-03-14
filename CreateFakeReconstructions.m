clear; clc;

load('xcat_reduced.mat');

slice = 5;
xcat_density_slice = xcat_density(:,:,slice);
xcat_id_slice = xcat_id(:,:,slice);

energy_vec = csvread('spectrum.txt');
max_energy = max(energy_vec);

materials_id = unique(xcat_id);
LoadGParams;
for ii=1:length(materials_id)
    Materials_names{ii} = gParamStruct.dict_mat_to_dens_names{find(gParamStruct.dict_mat_to_dens_xcat_id==materials_id(ii))};
    Materials_names{ii}(Materials_names{ii}==' ') = [];
end

% [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Full_Materials_names_clean);
[per_bin_linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names);
[per_energy_linear_atten_mat, ~] = PerEnergyMaterialsAttenuations(Materials_names, 1:length(energy_vec));

attenuation_slices = zeros([size(xcat_density_slice),6]);
attenuation_full_spectrum = zeros([size(xcat_density_slice)]);
for row = 1:size(xcat_density_slice,1)
    for col = 1:size(xcat_density_slice,2)
        attenuation_slices(row,col,:) = xcat_density_slice(row,col)*per_bin_linear_atten_mat(materials_id==xcat_id_slice(row,col),:);

        attenuation_full_spectrum(row,col) = xcat_density_slice(row,col)*sum(per_energy_linear_atten_mat(materials_id==xcat_id_slice(row,col),:).*energy_vec');
    end
end
%%
figure;
subplot(1,2,1);
imshow(attenuation_full_spectrum,[]);


for ii=1:6
    subplot(2,6,ii+3+3*(ii>3));
    imshow(attenuation_slices(:,:,ii),[]);
end

%%


noise_level = 0.05;

attenuation_full_spectrum_noised = attenuation_full_spectrum + noise_level*rand(size(attenuation_full_spectrum));
attenuation_slices_noised = attenuation_slices + noise_level*rand(size(attenuation_slices));

%%
figure;
subplot(1,2,1);
imshow(attenuation_full_spectrum_noised,[]);


for ii=1:6
    subplot(2,6,ii+3+3*(ii>3));
    imshow(attenuation_slices_noised(:,:,ii),[]);
end


%%

save('fake_reconstruction.mat','attenuation_full_spectrum','attenuation_full_spectrum_noised', ...
                               'attenuation_slices','attenuation_slices_noised', ...
                               'xcat_density_slice','xcat_id_slice');
                           
%%
clear; clc;
load('fake_reconstruction.mat')
LoadGParams;

%% Old Method:

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Getting Attenuation vector:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full_center_energy = 50;
% full_mu_water_e     = interp1((1e3)*gParamStruct.water_atten(:,1),0.99857*gParamStruct.water_atten(:,2),full_center_energy);
% full_recovery_ct = 1e3*((attenuation_full_spectrum_noised-full_mu_water_e)/(full_mu_water_e));
% full_recovery_density = interp1(gParamStruct.ct_to_dens(:,1),gParamStruct.ct_to_dens(:,2),full_recovery_ct);
% 
% materials_recovered_old = zeros(size(full_recovery_density));
% for ii=1:size(gParamStruct.dict_mat_to_dens_values,1)
%     materials_recovered_old = materials_recovered_old+ii*(full_recovery_density>=gParamStruct.dict_mat_to_dens_values(ii,1)).*(full_recovery_density<gParamStruct.dict_mat_to_dens_values(ii,2));
% end
% 
% old_index_to_id_dict = [23,6,0,2,0,1,0,18];
% %
% % % Only air-water-bone:
% % materials_recovered_old_awb = 23*ones(size(full_recovery_density));
% % materials_recovered_old_awb(full_recovery_density>=gParamStruct.dict_mat_to_dens_values(1,2)) = 1;
% % materials_recovered_old_awb(full_recovery_density>=gParamStruct.dict_mat_to_dens_values(8,1)) = 26;
% % 
% % 
% 


%% New Method:

noise_level = 0.035;

attenuation_full_spectrum_noised = attenuation_full_spectrum + noise_level*rand(size(attenuation_full_spectrum));
attenuation_slices_noised = attenuation_slices + noise_level*rand(size(attenuation_slices));


bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];

active_bins = [1,2,3,4,5,6];

% Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater'};
% Full_Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater','Blood','Redmarrow'};
Materials_names = {'Air','Water','Muscle','dryribwater','Redmarrow','Iodineblood','Adiposetissue'};
index_to_id_dict = [23,    1,      2           38             31         18              6];

% [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Full_Materials_names_clean);
[linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names);

energy_vec = csvread('spectrum.txt');
[per_energy_linear_atten_mat, ~] = PerEnergyMaterialsAttenuations(Materials_names, 1:length(energy_vec));
old_linear_atten = sum(per_energy_linear_atten_mat(:,:).*repmat(energy_vec',7,1),2);

active_linear_atten_2d_mat = linear_atten_mat(:,active_bins);
active_bin_weights = bin_weights(active_bins);

reconstructed_bins_vec = reshape(attenuation_slices_noised,[],6);

% materials_recovered_new = reshape(knnsearch(active_linear_atten_2d_mat,reconstructed_bins_vec,'Distance','cityblock'),size(attenuation_slices_noised(:,:,1)));
% materials_recovered_new = reshape(WeightedNearestNeighbour(active_linear_atten_2d_mat,active_bin_weights,reconstructed_bins_vec),size(attenuation_slices_noised(:,:,1)));

factor_per_different_label = 0.5;
materials_recovered_new = reshape(SmartWeightedNearestNeighbour(active_linear_atten_2d_mat,active_bin_weights,reconstructed_bins_vec,factor_per_different_label),size(attenuation_slices_noised(:,:,1)));
materials_recovered_old = reshape(SmartWeightedNearestNeighbour(old_linear_atten,1,attenuation_full_spectrum_noised(:),factor_per_different_label),size(attenuation_slices_noised(:,:,1)));

% Compare:
% xcat_id_recovered_old = gParamStruct.dict_mat_to_dens_xcat_id(materials_recovered_old);
xcat_id_recovered_old = index_to_id_dict(materials_recovered_old);
xcat_id_recovered_new = index_to_id_dict(materials_recovered_new);

error_old = sum(1-(xcat_id_recovered_old==xcat_id_slice),'all')/numel(xcat_id_slice);
error_new = sum(1-(xcat_id_recovered_new==xcat_id_slice),'all')/numel(xcat_id_slice);

fprintf('Old ID Errors: %d%%\n',round(error_old*100));
fprintf('New ID Errors: %d%%\n',round(error_new*100));
figure;
subplot(2,3,1);imagesc(xcat_id_recovered_new);
subplot(2,3,4);imagesc(xcat_id_recovered_old);
subplot(1,2,2);imagesc(xcat_id_slice)
% %%
% fractional_recovered_old = materials_db_mat(xcat_id_recovered_old(:),:);
% fractional_densities_recovered_old = fractional_recovered_old.*repmat(full_recovery_density(:),1,size(fractional_recovered_old,2));
% 
% fractional_recovered_new = materials_db_mat(xcat_id_recovered_new(:),:);
% fractional_densities_recovered_new = fractional_recovered_new.*repmat(full_recovery_density(:),1,size(fractional_recovered_new,2));
% 
% fractional_densities_xcat_gt = reshape(permute(xcat_gt.xcat_fractional_mass,[2,3,4,1]),[],27).*repmat(xcat_gt.xcat_density(:),1,size(fractional_recovered_new,2));
% 
% active_materials = 1:27; %[1,8,14,19];
% 
% density_error_old = sum((fractional_densities_recovered_old(:,active_materials)-fractional_densities_xcat_gt(:,active_materials)).^2,'all');
% density_error_new = sum((fractional_densities_recovered_new(:,active_materials)-fractional_densities_xcat_gt(:,active_materials)).^2,'all');
% 
% fprintf('Old Density Errors: %d\n',density_error_old);
% fprintf('New Density Errors: %d\n',density_error_new);
% fprintf('Density Errors new/old: %d%%\n',round(100*density_error_new/density_error_old));

%% New Method Comparing

bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];

active_bins = [1,2,3,4,5,6];

% Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater'};
% Full_Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater','Blood','Redmarrow'};
Materials_names = {'Air','Water','Muscle','dryribwater','Redmarrow','Iodineblood','Adiposetissue'};
index_to_id_dict = [23,    1,      2           38             31         18              6];

% [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Full_Materials_names_clean);
[linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names);

energy_vec = csvread('spectrum.txt');
[per_energy_linear_atten_mat, ~] = PerEnergyMaterialsAttenuations(Materials_names, 1:length(energy_vec));
old_linear_atten = sum(per_energy_linear_atten_mat(:,:).*repmat(energy_vec',7,1),2);

active_linear_atten_2d_mat = linear_atten_mat(:,active_bins);
active_bin_weights = bin_weights(active_bins);

noise_level_vec = 0:0.01:0.2;

num_runs_per_noise = 1;

error_old_vec = zeros(num_runs_per_noise,length(noise_level_vec));
error_new_vec = zeros(num_runs_per_noise,length(noise_level_vec));

old_recon_database = zeros(60,60,length(noise_level_vec));
new_recon_database = zeros(60,60,length(noise_level_vec));

for run_num = 1:num_runs_per_noise
    for jj=1:length(noise_level_vec)

        noise_level = noise_level_vec(jj);

        attenuation_full_spectrum_noised = attenuation_full_spectrum + noise_level*max(attenuation_full_spectrum,[],[1,2])*(rand(size(attenuation_full_spectrum))-0.5);
        attenuation_slices_noised = attenuation_slices + noise_level*repmat(max(attenuation_slices,[],[1,2]),60,60,1).*(rand(size(attenuation_slices))-0.5);

        reconstructed_bins_vec = reshape(attenuation_slices_noised,[],6);

        % materials_recovered_new = reshape(knnsearch(active_linear_atten_2d_mat,reconstructed_bins_vec,'Distance','cityblock'),size(attenuation_slices_noised(:,:,1)));
        % materials_recovered_new = reshape(WeightedNearestNeighbour(active_linear_atten_2d_mat,active_bin_weights,reconstructed_bins_vec),size(attenuation_slices_noised(:,:,1)));

        factor_per_different_label = 0.5;
        materials_recovered_new = reshape(SmartWeightedNearestNeighbour(active_linear_atten_2d_mat,active_bin_weights,reconstructed_bins_vec,factor_per_different_label),size(attenuation_slices_noised(:,:,1)));
        materials_recovered_old = reshape(SmartWeightedNearestNeighbour(old_linear_atten,1,attenuation_full_spectrum_noised(:),factor_per_different_label),size(attenuation_slices_noised(:,:,1)));

        % Compare:
        xcat_id_recovered_old = index_to_id_dict(materials_recovered_old);
        xcat_id_recovered_new = index_to_id_dict(materials_recovered_new);

        if run_num==1
            old_recon_database(:,:,jj) = xcat_id_recovered_old;
            new_recon_database(:,:,jj) = xcat_id_recovered_new;
        end
        
        error_old_vec(run_num,jj) = sum(1-(xcat_id_recovered_old==xcat_id_slice),'all')/numel(xcat_id_slice);
        error_new_vec(run_num,jj) = sum(1-(xcat_id_recovered_new==xcat_id_slice),'all')/numel(xcat_id_slice);
    end
end
figure;
subplot(2,1,1);
plot(mean(error_old_vec,1)); hold on;
plot(mean(error_new_vec,1));
legend('Old','New');
subplot(2,1,2);
plot(std(error_old_vec,1)); hold on;
plot(std(error_new_vec,1));
legend('Old','New');

figure;
subplot(2,3,1);imagesc(xcat_id_recovered_new);
subplot(2,3,4);imagesc(xcat_id_recovered_old);
subplot(1,2,2);imagesc(xcat_id_slice)
%%
noise_level = noise_level_vec(21);

attenuation_full_spectrum_noised = attenuation_full_spectrum + noise_level*(rand(size(attenuation_full_spectrum))-0.5);
attenuation_slices_noised = attenuation_slices + noise_level*(rand(size(attenuation_slices))-0.5);
figure;
subplot(1,2,1);
imshow(attenuation_full_spectrum_noised,[]);
title(sprintf('PSNR: %.2f',psnr(attenuation_full_spectrum_noised,attenuation_full_spectrum)));
for ii=1:6
    subplot(2,6,ii+3+3*(ii>3));
    imshow(attenuation_slices_noised(:,:,ii),[]);
    title(sprintf('PSNR: %.2f',psnr(attenuation_slices_noised(:,:,ii),attenuation_slices(:,:,ii))));
%     title(sprintf('Min: %.2f. Max: %.2f',min(attenuation_slices_noised(:,:,ii),[],'all'),max(attenuation_slices_noised(:,:,ii),[],'all')));
end


%%
%prepare figure and guidata struct
h=struct;
h.f=figure;
h.ax=axes('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.1 0.1 0.6 0.8]);
h.slider=uicontrol('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.8 0.1 0.1 0.8],...
    'Style','Slider',...
    'BackgroundColor',[1 1 1],...
    'Min',1,'Max',length(noise_level_vec),'Value',1,...
    'Callback',@sliderCallback);
%store image database to the guidata struct as well
% h.database=attenuation_full_spectrum;
h.database_gt  = xcat_id_slice;
h.database_old = old_recon_database;
h.database_new = new_recon_database;
guidata(h.f,h)
%trigger a callback
sliderCallback(h.slider)
function sliderCallback(hObject,~)
    h=guidata(hObject);
%     noise_level=(round(get(hObject,'Value'))-1)*0.001;
%     noised_image = h.database + noise_level*rand(size(h.database));
%     imshow(noised_image,[]);
    noise_level=round(get(hObject,'Value'));
    subplot(1,4,1); imagesc(h.database_old(:,:,noise_level)); title('old');
    subplot(1,4,2); imagesc(h.database_new(:,:,noise_level)); title('new');
    subplot(1,4,3); imagesc(h.database_gt); title('GT');
    title(['Noise Value: ',num2str(noise_level)]);
end

% %%
% 
% y = mean(error_old_vec,1);
% x = 1:numel(y);
% std_dev = std(error_old_vec,1);
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'g');
% hold on;
% plot(x, y, 'r', 'LineWidth', 2);
% 
% y = mean(error_new_vec,1);
% x = 1:numel(y);
% std_dev = std(error_new_vec,1);
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'g');
% hold on;
% plot(x, y, 'r', 'LineWidth', 2);