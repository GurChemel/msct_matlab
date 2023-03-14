clear;clc;

% load('reconstruction_atten_per_bin0.mat')
load('220113\reconstruction_atten_anti_grid_Enabled.mat')
% load('220208\reconstruction_atten_anti_grid_Enabled.mat')
load('clean_att_slices.mat');
% load('220212\reconstruction_atten_anti_grid_Disabled.mat')

slice_to_plot = round(size(reconstruct_bin_0,3)/2);
per_bin_image(:,:,1) = reconstruct_bin_0(:,:,slice_to_plot);
per_bin_image(:,:,2) = reconstruct_bin_1(:,:,slice_to_plot);
per_bin_image(:,:,3) = reconstruct_bin_2(:,:,slice_to_plot);
per_bin_image(:,:,4) = reconstruct_bin_3(:,:,slice_to_plot);
per_bin_image(:,:,5) = reconstruct_bin_4(:,:,slice_to_plot);
per_bin_image(:,:,6) = reconstruct_bin_5(:,:,slice_to_plot);

figure
for ii=1:6
    subplot(2,3,ii);
%     imagesc(per_bin_image(:,:,ii),[0,1.3]); colormap('bone');
    imagesc(per_bin_image(:,:,ii)); colormap('bone');
    axis off;
    title(sprintf('Bin %d',ii));
end
%%
att_vals = unique(reshape(attenuation_slices(:,:,1),1,[]));
pixels_per_val = zeros(1,length(att_vals));
for ii=1:length(att_vals)
    pixels_per_val(ii) = sum(attenuation_slices(:,:,1)==att_vals(ii),'all');
end
std_mat = zeros(6,length(att_vals));
avg_mat = zeros(6,length(att_vals));
att_mat = zeros(6,length(att_vals));
figure(6)
for ii=1:6

    rec_vec = reshape(per_bin_image(:,:,ii),1,[]);
    att_vec = reshape(attenuation_slices(:,:,ii),1,[]);
    att_vals = unique(att_vec);
    att_mat(ii,:) = att_vals;
    
    for att_val = att_vals
        std_mat(ii,att_val == att_vals) = std(rec_vec(att_vec==att_val));
        avg_mat(ii,att_val == att_vals) = mean(rec_vec(att_vec==att_val));
    end
    
    subplot(2,3,ii);
    scatter(rec_vec,att_vec); hold on;
    max_val = max([rec_vec,att_vec]);
    plot([0,max_val],[0,max_val],'r'); hold off;
    pbaspect([1,1,1])
    xlim([0,max_val]); ylim([0,max_val]);
    xlabel('Reconstructed Attenuation');
    ylabel('Ground Truth Attenuation');
    title(sprintf('Bin %d',ii));
%     for kk=1:length(att_vals)
%         text(max_val,att_vals(kk),sprintf('STD %.2f',std_mat(ii,kk)))
%     end
end
figure(7)
subplot(2,1,1);
plot(att_mat',std_mat','LineWidth',4); %hold on;
xlim([0,1])
% plot(repmat(att_vals',1,2)',repmat([0,0.15],7,1)',':k','LineWidth',1); hold off;
% for kk=1:length(att_vals)
%     text(att_vals(kk),0.13+0.01*mod(kk,2),sprintf('#Pixels=%d',pixels_per_val(kk)))
% end
xlabel('GT Attenuation [g/cm^2]'); ylabel('STD'); legend({'Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6'})
xticks(round(att_vals,3));
set(gca,'FontSize',19);
title('Standard Deviation Per Attenuation Per Bin')
subplot(2,1,2);
plot(att_mat,std_mat,'LineWidth',4); hold on;
%%
figure;
plot(avg_mat,'r'); hold on;
plot(att_mat,'k');
%%
per_bin_image_vec = reshape(per_bin_image,[],6);
per_bin_image_vec_norm = per_bin_image_vec./mean(per_bin_image_vec);
pca_res = pca(per_bin_image_vec_norm');
pca_res_im = reshape(pca_res,60,60,[]);
figure(8)
for ii=1:5
    subplot(2,3,ii);
    imagesc(pca_res_im(:,:,ii)); colormap('bone');
    axis off;
    title(sprintf('PCA %d',ii));
end

%%
num_cols = 4;

%%
selected_bins = [4,6];

% energy_centers = [31.3090   41.5364   50.3864   58.9392   69.0782   90.8517];
energy_centers = [30.5000   50.5000   60.0000   69.0000   81.5000  104.5000];

energy_low  = energy_centers(selected_bins(1));
energy_high = energy_centers(selected_bins(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Reconstructions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reconstruct_low_freq  = per_bin_image(:,:,selected_bins(1));
reconstruct_high_freq = per_bin_image(:,:,selected_bins(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Reconstructions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% subplot(2,3,1); 
% title_str = sprintf('MULTI Low Energy - %.2f [KeV]',energy_low);
% imagesc(reconstruct_low_freq); title(title_str);
% colormap('bone'); pbaspect([1,1,1]);
% write_im_for_pres(reconstruct_low_freq,title_str)
% subplot(2,3,4); 
% title_str = sprintf('MULTI Low Energy - %.2f [KeV]',energy_high);
% imagesc(reconstruct_high_freq); title(title_str);
% colormap('bone'); pbaspect([1,1,1]);
% write_im_for_pres(reconstruct_high_freq,title_str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Attenuation Matrix:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gParams;

% mu_water_e_low  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_low);
% mu_water_e_high = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_high);
mu_water_e_low  = interp1((1e3)*blood_atten(:,1),blood_atten(:,3),energy_low);
mu_water_e_high = interp1((1e3)*blood_atten(:,1),blood_atten(:,3),energy_high);
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
subplot(2,num_cols,1);
imagesc(norm_a1,[0,1]); t=title(sprintf('DECT from 2 bins\nWater Coefficient'));
colormap('bone'); pbaspect([1,1,1]); axis off;
% write_im_for_pres(norm_a1,t.String)
subplot(2,num_cols,1+num_cols);
imagesc(norm_a2,[0,1]); t=title(sprintf('DECT from 2 bins\nBone Coefficient'));
colormap('bone'); pbaspect([1,1,1]); axis off;
% write_im_for_pres(norm_a2,t.String)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New Method:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mu_water, ~] = PerBinMaterialsAttenuations({'Water'});
[mu_bone, ~] = PerBinMaterialsAttenuations({'dryribwater'});
% mu_water = interp1((1e3)*blood_atten(:,1),blood_atten(:,3),energy_centers)';
% mu_bone = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_centers)';
mu_matrix_full = [mu_water,mu_bone];
% weight_vector = [0,0,0.5,0.6,0.4,0.5];
% weight_vector = [0,0,0,.56,0,0.6];
% weight_vector = (1./mu_water).^0.7;
% weight_vector = (1./std(reshape(per_bin_image,[],6)));
weight_vector = ones(1,6);
% rcond(mu_matrix_full'*diag(weight_vector)*mu_matrix_full)
% inv_mu_matrix_full = ((mu_matrix_full'*diag(weight_vector)*mu_matrix_full)^-1)*mu_matrix_full'*diag(weight_vector);
inv_mu_matrix_full = ((mu_matrix_full'*diag(weight_vector)*mu_matrix_full)^-1)*(mu_matrix_full)';

a_mat = inv_mu_matrix_full*(reshape(per_bin_image,[],6)');
a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,num_cols,2);
imagesc(norm_a1,[0,1]); t=title(sprintf('Multispectal No Weight Vector\nWater Coefficient'));
colormap('bone'); pbaspect([1,1,1]); axis off;
subplot(2,num_cols,2+num_cols);
imagesc(norm_a2,[0,1]); t=title(sprintf('Multispectal No Weight Vector\nBone Coefficient'));
colormap('bone'); pbaspect([1,1,1]); axis off;
%%
% weight_vector = [0,0,0.5,0.6,0.4,0.5];
% weight_vector = [0,0,0,.56,0,0.6];
% weight_vector = (1./mu_water).^0.7;
% weight_vector = (1./std(reshape(per_bin_image,[],6)));
weight_vector = (1./(std_mat(:,2)'));

% rcond(mu_matrix_full'*diag(weight_vector)*mu_matrix_full)
% inv_mu_matrix_full = ((mu_matrix_full'*diag(weight_vector)*mu_matrix_full)^-1)*(diag(weight_vector)*mu_matrix_full)';
inv_mu_matrix_full = ((mu_matrix_full'*diag(weight_vector)*mu_matrix_full)^-1)*(mu_matrix_full)';

a_mat = inv_mu_matrix_full*(reshape(per_bin_image,[],6)');
a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,num_cols,3);
imagesc(norm_a1,[0,1]); t=title(sprintf('Multispectal (STD)^{-1} Weight Vector\nWater Coefficient'));
colormap('bone'); pbaspect([1,1,1]); axis off;
subplot(2,num_cols,3+num_cols);
imagesc(norm_a2,[0,1]); t=title(sprintf('Multispectal (STD)^{-1} Weight Vector\nBone Coefficient'));
colormap('bone'); pbaspect([1,1,1]); axis off;

%%
% weight_vector = [0,0,0.5,0.6,0.4,0.5];
% weight_vector = [0,0,0,0.5,0.2,0.1];
% weight_vector = (1./mu_water).^0.7;
% weight_vector = (1./std(reshape(per_bin_image,[],6))).^1.3;
weight_vector = (1./((std_mat(:,2)').^2));

% rcond(mu_matrix_full'*diag(weight_vector)*mu_matrix_full)
% inv_mu_matrix_full = ((mu_matrix_full'*diag(weight_vector)*mu_matrix_full)^-1)*(diag(weight_vector)*mu_matrix_full)';
inv_mu_matrix_full = ((mu_matrix_full'*diag(weight_vector)*mu_matrix_full)^-1)*(mu_matrix_full)';

a_mat = inv_mu_matrix_full*(reshape(per_bin_image,[],6)');
a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,num_cols,4);
imagesc(norm_a1,[0,1]); t=title(sprintf('Multispectal (STD^{2})^{-1} Weight Vector\nWater Coefficient'));
colormap('bone'); pbaspect([1,1,1]); axis off;
write_im_for_pres(norm_a1,'MECT_STD_pow_2_water')
subplot(2,num_cols,4+num_cols);
imagesc(norm_a2,[0,1]); t=title(sprintf('Multispectal (STD^{2})^{-1} Weight Vector\nBone Coefficient'));
colormap('bone'); pbaspect([1,1,1]); axis off;
write_im_for_pres(norm_a2,'MECT_STD_pow_2_bone')

% %%
% weight_vector = [0.2,0.5,0.3,0.8,0.8,0.8]*20;
% % weight_vector = [0,0,0,.56,0,0.6];
% % weight_vector = (1./mu_water).^0.7;
% % weight_vector = (1./std(reshape(per_bin_image,[],6))).^1.3;
% % rcond(mu_matrix_full'*diag(weight_vector)*mu_matrix_full)
% inv_mu_matrix_full = (((mu_matrix_full)'*diag(weight_vector)*mu_matrix_full)^-1)*(diag(weight_vector)*mu_matrix_full)';
% 
% a_mat = inv_mu_matrix_full*(reshape(per_bin_image,[],6)');
% a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
% a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));
% 
% norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
% norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plotting Decompositions:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,num_cols,5);
% imagesc(norm_a1,[0,1]); t=title(sprintf('Multispectal Custom Weight Vector\nWater Coefficient'));
% colormap('bone'); pbaspect([1,1,1]); axis off;
% subplot(2,num_cols,5+num_cols);
% imagesc(norm_a2,[0,1]); t=title(sprintf('Multispectal Custom Weight Vector\nBone Coefficient'));
% colormap('bone'); pbaspect([1,1,1]); axis off;
