clear; clc;

load('220201\xcat_reduced.mat');

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
%         attenuation_slices(row,col,:) = xcat_density_slice(row,col)*per_bin_linear_atten_mat(materials_id==xcat_id_slice(row,col),:);
%         attenuation_full_spectrum(row,col) = xcat_density_slice(row,col)*sum(per_energy_linear_atten_mat(materials_id==xcat_id_slice(row,col),:).*energy_vec');
        attenuation_slices(row,col,:) = per_bin_linear_atten_mat(materials_id==xcat_id_slice(row,col),:);
        attenuation_full_spectrum(row,col) = sum(per_energy_linear_atten_mat(materials_id==xcat_id_slice(row,col),:).*energy_vec');
    end
end

%%
clc;
load('selected_std.mat');

% load('220201\reconstruction_atten_anti_grid_Enabled.mat');
% attenuation_slices_noised(:,:,1) = reconstruct_bin_0(:,:,slice);
% attenuation_slices_noised(:,:,2) = reconstruct_bin_1(:,:,slice);
% attenuation_slices_noised(:,:,3) = reconstruct_bin_2(:,:,slice);
% attenuation_slices_noised(:,:,4) = reconstruct_bin_3(:,:,slice);
% attenuation_slices_noised(:,:,5) = reconstruct_bin_4(:,:,slice);
% attenuation_slices_noised(:,:,6) = reconstruct_bin_5(:,:,slice);

noise_level = 0.25;
noise_factor = permute([1,1.1,1.35,1.5,1.8,2],[1,3,2]);
attenuation_slices_noised = attenuation_slices + noise_level*repmat(noise_factor.*max(attenuation_slices,[],[1,2]),60,60,1).*(rand(size(attenuation_slices))-0.5);

fake_std_vec = std(reshape(attenuation_slices_noised(repmat(xcat_id_slice,1,1,6)==6),[],6));
fprintf('Simulation Tissue STD Vec = ');fprintf('%.4f, ',simulation_tissue_std_per_bin);fprintf('\b\b.\n');
fprintf('Fake Image Tissue STD Vec = ');fprintf('%.4f, ',fake_std_vec);fprintf('\b\b.\n');

%% Recon Method Comparing

addpath('GCMex/')

LoadGParams;
bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];

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

% active_linear_atten_2d_mat = linear_atten_mat(:,active_bins);
% active_bin_weights = bin_weights(active_bins);

% noise_level_vec = 0.2:0.1:0.4;
% noise_level_vec = 0:0.01:0.2;
noise_level_vec = 0:0.4:0.8;

num_runs_per_noise = 10;

error_old_vec = zeros(num_runs_per_noise,length(noise_level_vec));
error_mid_vec = zeros(num_runs_per_noise,length(noise_level_vec));
error_new_vec = zeros(num_runs_per_noise,length(noise_level_vec));
error_gct_vec = zeros(num_runs_per_noise,length(noise_level_vec));

psnr_vec = zeros(num_runs_per_noise,length(noise_level_vec));

old_recon_database = zeros(60,60,length(noise_level_vec));
mid_recon_database = zeros(60,60,length(noise_level_vec));
new_recon_database = zeros(60,60,length(noise_level_vec));
gct_recon_database = zeros(60,60,length(noise_level_vec));

for run_num = 1:num_runs_per_noise
    for jj=1:length(noise_level_vec)

        noise_level = noise_level_vec(jj);

        attenuation_full_spectrum_noised = attenuation_full_spectrum + noise_level*max(attenuation_full_spectrum,[],[1,2])*(rand(size(attenuation_full_spectrum))-0.5);
        attenuation_slices_noised = attenuation_slices + noise_level*repmat(max(attenuation_slices,[],[1,2]),60,60,1).*(rand(size(attenuation_slices))-0.5);

        reconstructed_bins_vec = reshape(attenuation_slices_noised,[],6);

        % materials_recovered_new = reshape(knnsearch(active_linear_atten_2d_mat,reconstructed_bins_vec,'Distance','cityblock'),size(attenuation_slices_noised(:,:,1)));
        % materials_recovered_new = reshape(WeightedNearestNeighbour(active_linear_atten_2d_mat,active_bin_weights,reconstructed_bins_vec),size(attenuation_slices_noised(:,:,1)));

        factor_per_different_label = 0.5;
        smothness_factor = 0.03;
%         materials_recovered_old = reshape(SmartWeightedNearestNeighbour(full_spectrum_linear_atten,1,attenuation_full_spectrum_noised(:),factor_per_different_label),size(attenuation_slices_noised(:,:,1)));
        materials_recovered_old = GraphCutAccumWrapper(attenuation_full_spectrum_noised, smothness_factor);
        materials_recovered_mid = reshape(SmartWeightedNearestNeighbour(linear_atten_mat(:,[4,6]),bin_weights([4,6]),reconstructed_bins_vec(:,[4,6]),factor_per_different_label),size(attenuation_slices_noised(:,:,1)));
        materials_recovered_new = reshape(SmartWeightedNearestNeighbour(linear_atten_mat,bin_weights,reconstructed_bins_vec,factor_per_different_label),size(attenuation_slices_noised(:,:,1)));
        materials_recovered_gct = GraphCutPerBinWrapper(attenuation_slices_noised, smothness_factor);
        
        % Compare:
        xcat_id_recovered_old = index_to_id_dict(materials_recovered_old);
        xcat_id_recovered_mid = index_to_id_dict(materials_recovered_mid);
        xcat_id_recovered_new = index_to_id_dict(materials_recovered_new);
        xcat_id_recovered_gct = index_to_id_dict(materials_recovered_gct);

        if run_num==1
            old_recon_database(:,:,jj) = xcat_id_recovered_old;
            mid_recon_database(:,:,jj) = xcat_id_recovered_mid;
            new_recon_database(:,:,jj) = xcat_id_recovered_new;
            gct_recon_database(:,:,jj) = xcat_id_recovered_gct;
        end
        
        psnr_vec(run_num,jj) = psnr(attenuation_full_spectrum_noised,attenuation_full_spectrum);
        
        error_old_vec(run_num,jj) = sum(1-(xcat_id_recovered_old==xcat_id_slice),'all')/numel(xcat_id_slice);
        error_mid_vec(run_num,jj) = sum(1-(xcat_id_recovered_mid==xcat_id_slice),'all')/numel(xcat_id_slice);
        error_new_vec(run_num,jj) = sum(1-(xcat_id_recovered_new==xcat_id_slice),'all')/numel(xcat_id_slice);
        error_gct_vec(run_num,jj) = sum(1-(xcat_id_recovered_gct==xcat_id_slice),'all')/numel(xcat_id_slice);
        
        if run_num == 1
            if jj == 1
                figure;
            end
            subplot(length(noise_level_vec),4,4*(jj-1)+1); imagesc(xcat_id_recovered_old); pbaspect([1,1,1]); title('Full Spectrum Detectors');
            subplot(length(noise_level_vec),4,4*(jj-1)+2); imagesc(xcat_id_recovered_new); pbaspect([1,1,1]); title('Multi Spectral Detectors');
            subplot(length(noise_level_vec),4,4*(jj-1)+3); imagesc(xcat_id_recovered_gct); pbaspect([1,1,1]); title('Multi Spectral Detectors - GraphCut');
            subplot(length(noise_level_vec),4,4*(jj-1)+4); imagesc(xcat_id_slice); pbaspect([1,1,1]); title('Ground Truth');
        end

    end
end
%%
psnr_vec = mean(psnr_vec,1);

figure;
subplot(2,1,1);
plot(psnr_vec,mean(error_old_vec,1)); hold on;
plot(psnr_vec,mean(error_mid_vec,1));
plot(psnr_vec,mean(error_new_vec,1));
plot(psnr_vec,mean(error_gct_vec,1));
legend('Old','Mid','New','GraphCut');
subplot(2,1,2);
plot(psnr_vec,std(error_old_vec,1)); hold on;
plot(psnr_vec,std(error_mid_vec,1));
plot(psnr_vec,std(error_new_vec,1));
plot(psnr_vec,std(error_gct_vec,1));
legend('Old','Mid','New','GraphCut');

% figure;
% subplot(2,2,1);imagesc(xcat_id_recovered_old);
% subplot(2,2,2);imagesc(xcat_id_recovered_mid);
% subplot(2,2,3);imagesc(xcat_id_recovered_new);
% subplot(2,2,4);imagesc(xcat_id_slice)
%
% 
% round = 1;
% 
% if round == 0
%     figure;
% end
% subplot(2,4,1+4*round);imagesc(xcat_id_recovered_old); title('Full Spectrum Detectors');
% subplot(2,4,2+4*round);imagesc(xcat_id_recovered_mid); title('2 Energy Bins from 6');
% subplot(2,4,3+4*round);imagesc(xcat_id_recovered_new); title('6 Energy Bins from 6');
% subplot(2,4,4+4*round);imagesc(xcat_id_slice); title('Ground Truth');
%%
figure;
subplot(1,4,1); imagesc(xcat_id_recovered_old); pbaspect([1,1,1]); title('Full Spectrum Detectors');
subplot(1,4,2); imagesc(xcat_id_recovered_new); pbaspect([1,1,1]); title('Multi Spectral Detectors');
subplot(1,4,3); imagesc(xcat_id_recovered_gct); pbaspect([1,1,1]); title('Multi Spectral Detectors - GraphCut');
subplot(1,4,4); imagesc(xcat_id_slice); pbaspect([1,1,1]); title('Ground Truth');

%%
figure;
tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
nexttile; imagesc(xcat_id_recovered_old); pbaspect([1,1,1]); xlabel('(1)'); set(gca,'xticklabel',[],'yticklabel',[],'FontSize',16);
nexttile; imagesc(xcat_id_recovered_new); pbaspect([1,1,1]); xlabel('(2)'); set(gca,'xticklabel',[],'yticklabel',[],'FontSize',16);
nexttile; imagesc(xcat_id_recovered_gct); pbaspect([1,1,1]); xlabel('(3)'); set(gca,'xticklabel',[],'yticklabel',[],'FontSize',16);
nexttile; imagesc(xcat_id_slice);         pbaspect([1,1,1]); xlabel('(4)'); set(gca,'xticklabel',[],'yticklabel',[],'FontSize',16);

%%
figure;
% subplot(1,2,1);
t1 = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile(t1);
imshow(attenuation_full_spectrum_noised,[]);
xlabel('(0)'); set(gca,'xticklabel',[],'yticklabel',[],'FontSize',16);
t2=tiledlayout(t1,2,3,'TileSpacing','Compact','Padding','Compact');
for ii=1:6
    t2.Layout.Tile = 2;
    t2.Layout.TileSpan = [1 1];
    nexttile(t2);
%     subplot(2,6,ii+3+3*(ii>3));
   
    im_for_saving = rescale(attenuation_slices_noised(:,:,ii));
    imwrite(kron(im_for_saving,ones(3)),sprintf('ImagesForPresentations/labeling_noised_%d.png',ii))
    imshow(im_for_saving);
    xlabel(sprintf('(%d)',ii)); set(gca,'xticklabel',[],'yticklabel',[],'FontSize',16);
end

%%
% noise_level = noise_level_vec(2);

attenuation_full_spectrum_noised = attenuation_full_spectrum + noise_level*(rand(size(attenuation_full_spectrum))-0.5);
attenuation_slices_noised = attenuation_slices + noise_level*(rand(size(attenuation_slices))-0.5);
figure;
subplot(1,2,1);
imshow(attenuation_full_spectrum_noised,[]);
title(sprintf('Full Spectrum Attenuation Map. PSNR: %.2f',psnr(attenuation_full_spectrum_noised,attenuation_full_spectrum)));
for ii=1:6
    subplot(2,6,ii+3+3*(ii>3));
    imshow(attenuation_slices_noised(:,:,ii),[]);
    title(sprintf('Bin %d. PSNR: %.2f',ii, psnr(attenuation_slices_noised(:,:,ii),attenuation_slices(:,:,ii))));
%     title(sprintf('Min: %.2f. Max: %.2f',min(attenuation_slices_noised(:,:,ii),[],'all'),max(attenuation_slices_noised(:,:,ii),[],'all')));
end


%%
% %prepare figure and guidata struct
% h=struct;
% h.f=figure;
% h.ax=axes('Parent',h.f,...
%     'Units','Normalized',...
%     'Position',[0.1 0.1 0.6 0.8]);
% h.slider=uicontrol('Parent',h.f,...
%     'Units','Normalized',...
%     'Position',[0.8 0.1 0.1 0.8],...
%     'Style','Slider',...
%     'BackgroundColor',[1 1 1],...
%     'Min',1,'Max',length(noise_level_vec),'Value',1,...
%     'Callback',@sliderCallback);
% %store image database to the guidata struct as well
% % h.database=attenuation_full_spectrum;
% h.database_gt  = xcat_id_slice;
% h.database_old = old_recon_database;
% h.database_new = new_recon_database;
% guidata(h.f,h)
% %trigger a callback
% sliderCallback(h.slider)
% function sliderCallback(hObject,~)
%     h=guidata(hObject);
% %     noise_level=(round(get(hObject,'Value'))-1)*0.001;
% %     noised_image = h.database + noise_level*rand(size(h.database));
% %     imshow(noised_image,[]);
%     noise_level=round(get(hObject,'Value'));
%     subplot(1,4,1); imagesc(h.database_old(:,:,noise_level)); title('old');
%     subplot(1,4,2); imagesc(h.database_new(:,:,noise_level)); title('new');
%     subplot(1,4,3); imagesc(h.database_gt); title('GT');
%     title(['Noise Value: ',num2str(noise_level)]);
% end

%% Recon Method Comparing

LoadGParams;
bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];

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

% active_linear_atten_2d_mat = linear_atten_mat(:,active_bins);
% active_bin_weights = bin_weights(active_bins);

noise_level_vec = noise_level_vec(end); %0.5;
smothness_factor_vec = 0.01; %0.00:0.01:0.04;
num_runs_per_noise = 10;

error_old_vec = zeros(num_runs_per_noise,length(noise_level_vec));
error_mid_vec = zeros(num_runs_per_noise,length(noise_level_vec));
error_new_vec = zeros(num_runs_per_noise,length(noise_level_vec));
error_gct_vec = zeros(num_runs_per_noise,length(noise_level_vec));

psnr_vec = zeros(num_runs_per_noise,length(noise_level_vec));

old_recon_database = zeros(60,60,length(noise_level_vec));
mid_recon_database = zeros(60,60,length(noise_level_vec));
new_recon_database = zeros(60,60,length(noise_level_vec));
gct_recon_database = zeros(60,60,length(noise_level_vec));

noise_level = 0.5;
attenuation_full_spectrum_noised = attenuation_full_spectrum + noise_level*max(attenuation_full_spectrum,[],[1,2])*(rand(size(attenuation_full_spectrum))-0.5);
attenuation_slices_noised = attenuation_slices + noise_level*repmat(max(attenuation_slices,[],[1,2]),60,60,1).*(rand(size(attenuation_slices))-0.5);

for run_num = 1:num_runs_per_noise
    for jj=1:length(smothness_factor_vec)

        factor_per_different_label = 0.5;
        smothness_factor = smothness_factor_vec(jj);
        materials_recovered_old = GraphCutAccumWrapper(attenuation_full_spectrum_noised, smothness_factor);
        materials_recovered_gct = GraphCutPerBinWrapper(attenuation_slices_noised, smothness_factor);
        
        % Compare:
        xcat_id_recovered_old = index_to_id_dict(materials_recovered_old);
        xcat_id_recovered_gct = index_to_id_dict(materials_recovered_gct);

        if run_num == 1
            if jj == 1
                figure;
            end
            subplot(2,length(smothness_factor_vec),jj); imagesc(xcat_id_recovered_old); pbaspect([1,1,1]); title('Full Spectrum Detectors');
            subplot(2,length(smothness_factor_vec),jj+length(smothness_factor_vec)); imagesc(xcat_id_recovered_gct); pbaspect([1,1,1]); title('Multi Spectral Detectors - GraphCut');
        end

    end
end
%%

imwrite( ind2rgb(im2uint8(mat2gray(kron(xcat_id_slice,ones(3)))), parula(256)), 'ImagesForPresentations/labeling_id_gt.png')
imwrite( ind2rgb(im2uint8(mat2gray(kron(xcat_id_recovered_gct,ones(3)))), parula(256)), 'ImagesForPresentations/labeling_id_gct.png')