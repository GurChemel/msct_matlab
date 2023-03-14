clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading:
load('xcat_reduced.mat','xcat_density');
load('reconstruction_full_spectrum.mat','reconstruct_with_asg');
load('reconstruction_atten_per_bin0.mat');
load('init_density_compare.mat')

slice = 5;
gt_density_image     = xcat_density(:,:,slice);
full_spectrum_image  = reconstruct_with_asg(:,:,slice);
per_bin_image(:,:,1) = reconstruct_bin_0(:,:,slice);
per_bin_image(:,:,2) = reconstruct_bin_1(:,:,slice);
per_bin_image(:,:,3) = reconstruct_bin_2(:,:,slice);
per_bin_image(:,:,4) = reconstruct_bin_3(:,:,slice);
per_bin_image(:,:,5) = reconstruct_bin_4(:,:,slice);
per_bin_image(:,:,6) = reconstruct_bin_5(:,:,slice);

% energy_centers = [31.3090   41.5364   50.3864   58.9392   69.0782   90.8517].';
energy_centers = [31.3   41.536   50.386   58.939   69.078   90.852].';
full_center_energy = 50;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plotting Reconstructions:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% subplot(1,3,1)
% imagesc(full_spectrum_image);
% title(sprintf('Recontruction With Full Spectrum'));
% colormap('bone'); pbaspect([1,1,1]);
% for ii=1:6
%     subplot(2,6,ii+3+(ii>3)*3); 
%     imagesc(squeeze(per_bin_image(:,:,ii))); title(sprintf('Recontruction With Energy - %.2f [KeV]',energy_centers(ii)));
%     colormap('bone'); pbaspect([1,1,1]);
% end

% figure;
% imagesc(reshape(per_bin_image,60,[]))
% colormap('bone'); pbaspect([6,1,1]);
% title('Same scale per-bin attenuations');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Attenuation vector:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gParams;
full_mu_water_e     = interp1((1e3)*water_atten(:,1),0.99857*water_atten(:,2),full_center_energy);
per_bin_mu_water_e  = interp1((1e3)*water_atten(:,1),0.99857*water_atten(:,2),energy_centers);

full_recovery_ct = 1e3*((full_spectrum_image-full_mu_water_e)/(full_mu_water_e));
full_recovery_density = interp1(ct_to_dens(:,1),ct_to_dens(:,2),full_recovery_ct);

per_bin_recovery_density = zeros(size(per_bin_image));
for ii=1:6
    per_bin_recovery_ct = 1e3*((per_bin_image(:,:,ii)-per_bin_mu_water_e(ii))/(per_bin_mu_water_e(ii)));
    per_bin_recovery_density(:,:,ii) = interp1(ct_to_dens(:,1),ct_to_dens(:,2),per_bin_recovery_ct);
end

% per_bin_recovery_density = squeeze(density_per_bin(:,:,slice,:));
figure;
for ii=1:6
    subplot(2,6,ii)
    imagesc(per_bin_recovery_density(:,:,ii),[0 3.5])
    title(sprintf('Matlab Recovery. Bin: %d',ii));
    pbaspect([1,1,1]); colorbar();
    subplot(2,6,6+ii)
    imagesc(density_per_bin(:,:,slice,ii),[0 3.5])
    title(sprintf('Python Recovery. Bin: %d',ii));
    pbaspect([1,1,1]); colorbar();
    
    % Taking Python Values:
    per_bin_recovery_density(:,:,ii) = density_per_bin(:,:,slice,ii);
end

% tmp =  (per_bin_recovery_density-squeeze(density_per_bin(:,:,slice,:)));
% figure;
% for ii=1:6
%     subplot(1,6,ii)
%     imagesc(tmp(:,:,ii));
%     pbaspect([1,1,1]); colorbar();
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Density Reconstructions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% subplot(1,3,1)
% imagesc(full_recovery_density);
% title(sprintf('Recontruction With Full Spectrum'));
% pbaspect([1,1,1]);
% for ii=1:6
%     subplot(2,6,ii+3+(ii>3)*3); 
%     imagesc(squeeze(per_bin_recovery_density(:,:,ii))); title(sprintf('Recontruction With Energy - %.2f [KeV]',energy_centers(ii)));
%     pbaspect([1,1,1]);
% end
% figure;
% imagesc(reshape(per_bin_recovery_density,60,[]))
% pbaspect([6,1,1]);
% title('Same scale per-bin attenuations');

%%
% energy_vec = csvread('spectrum.txt')';
% full_water_atten = interp1((1e3)*water_atten(:,1),water_atten(:,3),1:150);
% water_dec_energy_vec = energy_vec.*exp(-full_water_atten);
% 
% ratio_vec = water_dec_energy_vec(round(energy_centers))./energy_vec(round(energy_centers));

% combined_per_bin_recovered_density = sum(per_bin_recovery_density,3)/6;

% ratio_vec = 1./[761828, 798509, 818354, 895294, 917473, 1009081];
% 
% % covariance_matrix = diag(var(reshape(per_bin_recovery_density,[],6)));
% covariance_matrix = diag(sqrt(ratio_vec));
% multiplexing_matrix = ones(6,1);
% wighted_least_squares_mult_matrix = (multiplexing_matrix.'*(covariance_matrix)^-1*multiplexing_matrix)\(multiplexing_matrix.')*((covariance_matrix)^-1);
% combined_per_bin_recovered_density = reshape(wighted_least_squares_mult_matrix*(reshape(per_bin_recovery_density,[],6).'),60,60);

figure;
for ii=1:6
    subplot(1,6,ii)
    imagesc(per_bin_recovery_density(:,:,ii));
    pbaspect([1,1,1]); colorbar();
end
combined_per_bin_recovered_density = reshape(median(reshape(per_bin_recovery_density,[],6).'),60,60);


clim_vec = [0,1.5];
figure;
subplot(1,3,1)
imagesc(full_recovery_density);
title(sprintf('Recontruction With Full Spectrum'));
pbaspect([1,1,1]); colorbar();
subplot(1,3,2)
imagesc(combined_per_bin_recovered_density)
title(sprintf('Recontruction With Per-Bin Summation'));
pbaspect([1,1,1]); colorbar();
subplot(1,3,3)
imagesc(gt_density_image)
title(sprintf('Ground truth'));
pbaspect([1,1,1]); colorbar();

%%
% %%
% figure;
% for ii=1:6
%     subplot(2,6,ii+3+(ii>3)*3); 
%     imagesc(squeeze(per_bin_image(:,:,ii))); title(sprintf('Recontruction With Energy - %.2f [KeV]',energy_centers(ii)));
%     colormap('bone'); pbaspect([1,1,1]);
% end


%%
% % close all;
% 
% air_thresh = 0.28;
% 
% mul_init_no_air = mul_init_density;
% mul_init_no_air(mul_init_no_air<air_thresh) = 0;
% 
% clim_vec = [0,1.5];
% figure;
% subplot(2,3,1)
% imagesc(lin_init_density(:,:,slice),clim_vec);
% title(sprintf('Full Spectrum. MSE = %.4f. PSNR = %.2f',immse(lin_init_density(:,:,slice),gt_density(:,:,slice)),psnr(lin_init_density(:,:,slice),gt_density(:,:,slice))));
% pbaspect([1,1,1]); colorbar();
% subplot(2,3,2)
% imagesc(mul_init_no_air(:,:,slice),clim_vec)
% title(sprintf('Per-Bin WLS. MSE = %.4f. PSNR = %.2f',immse(mul_init_no_air(:,:,slice),gt_density(:,:,slice)),psnr(mul_init_no_air(:,:,slice),gt_density(:,:,slice))));
% pbaspect([1,1,1]); colorbar();
% subplot(2,3,3)
% imagesc(gt_density(:,:,slice),clim_vec)
% title(sprintf('Ground truth'));
% pbaspect([1,1,1]); colorbar();
% 
% subplot(2,3,4)
% imagesc(lin_init_density(:,:,slice));
% title(sprintf('Recontruction With Full Spectrum'));
% pbaspect([1,1,1]); colorbar();
% subplot(2,3,5)
% imagesc(mul_init_no_air(:,:,slice))
% title(sprintf('Recontruction With Per-Bin Summation'));
% pbaspect([1,1,1]); colorbar();
% subplot(2,3,6)
% imagesc(gt_density(:,:,slice))
% title(sprintf('Ground truth'));
% pbaspect([1,1,1]); colorbar();

%%
% close all;

air_thresh = 0;

% % complicated_ratio_vec = ratio_vec.*linspace(20,1,6);
% % complicated_ratio_vec = sqrt(1./[761828, 798509, 818354, 895294, 917473, 1009081]);
% complicated_ratio_vec = sqrt(linspace(20,1,6));
% covariance_matrix = diag((complicated_ratio_vec));
% multiplexing_matrix = ones(6,1);
% wighted_least_squares_mult_matrix = (multiplexing_matrix.'*(covariance_matrix)^-1*multiplexing_matrix)\(multiplexing_matrix.')*((covariance_matrix)^-1);
% combined_per_bin_recovered_density = reshape(wighted_least_squares_mult_matrix*(reshape(per_bin_recovery_density,[],6).'),60,60);

combined_per_bin_recovered_density_no_air = combined_per_bin_recovered_density;
combined_per_bin_recovered_density_no_air(combined_per_bin_recovered_density_no_air<air_thresh) = 0;

clim_vec = [0,1.5];
figure;
subplot(2,3,1)
imagesc(lin_init_density(:,:,slice),clim_vec);
title(sprintf('Full Spectrum. MSE = %.4f. PSNR = %.2f',immse(lin_init_density(:,:,slice),gt_density(:,:,slice)),psnr(lin_init_density(:,:,slice),gt_density(:,:,slice))));
pbaspect([1,1,1]); colorbar();
subplot(2,3,2)
imagesc(combined_per_bin_recovered_density_no_air,clim_vec)
title(sprintf('Per-Bin WLS. MSE = %.4f. PSNR = %.2f',immse(combined_per_bin_recovered_density_no_air,gt_density(:,:,slice)),psnr(combined_per_bin_recovered_density_no_air,gt_density(:,:,slice))));
pbaspect([1,1,1]); colorbar();
subplot(2,3,3)
imagesc(gt_density(:,:,slice),clim_vec)
title(sprintf('Ground truth'));
pbaspect([1,1,1]); colorbar();

subplot(2,3,4)
imagesc(lin_init_density(:,:,slice));
title(sprintf('Recontruction With Full Spectrum'));
pbaspect([1,1,1]); colorbar();
subplot(2,3,5)
imagesc(combined_per_bin_recovered_density_no_air)
title(sprintf('Recontruction With Per-Bin Summation'));
pbaspect([1,1,1]); colorbar();
subplot(2,3,6)
imagesc(gt_density(:,:,slice))
title(sprintf('Ground truth'));
pbaspect([1,1,1]); colorbar();

%%
% 
% mat_id = zeros(size(combined_per_bin_recovered_density));
% for xx=1:size(combined_per_bin_recovered_density,1)
%     for yy=1:size(combined_per_bin_recovered_density,2)
%         mat_id(xx,yy) = find((combined_per_bin_recovered_density(xx,yy)>=dict_mat_to_dens_values(:,1)) .* (combined_per_bin_recovered_density(xx,yy)<=dict_mat_to_dens_values(:,2)),1);
%     end
% end
% 
% 
% slice_id = mat_id;
% uniqe_id_value = unique(slice_id);
% num_uniques = length(uniqe_id_value);
% simplified_id_map = slice_id;
% for ii=1:num_uniques
%     simplified_id_map(slice_id==uniqe_id_value(ii))=ii;
%     material_name = dict_mat_to_dens_names{uniqe_id_value(ii)};
%     material_name(material_name=='_')=' ';
%     simplified_legend{ii} = material_name;
% end
% 
% figure;
% cmap = hsv(num_uniques);
% image(simplified_id_map);
% colormap(cmap);
% hold on;
% for K = 1 : num_uniques; hidden_h(K) = surf(uint8(K-[1 1;1 1]), 'edgecolor', 'none'); end
% hold off
% uistack(hidden_h, 'bottom');
% legend(hidden_h, simplified_legend )
