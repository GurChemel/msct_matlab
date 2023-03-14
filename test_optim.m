clear; clc;
addpath('PhotonAttenuation\')

% xcat_db = load('xcat_reduced.mat');
% init_db = load('Optimization\Single_0_0init_vol.mat');
% final_0_db = load('Optimization\Single_0final_vol.mat');
% final_1_db = load('Optimization\Single_1final_vol.mat');

folder = '210819\';

xcat_db = load([folder,'xcat_reduced.mat']);
init_db = load([folder,'Single_0_0init_vol.mat']);
final_0_db = load([folder,'Single_0final_vol.mat']);
final_1_db = load([folder,'Single_1final_vol.mat']);

slice = 5;

figure;
subplot(2,3,1);
imagesc(init_db.init_vol_density(:,:,slice));
title(sprintf('Init Vol Density'));
pbaspect([1,1,1]); colorbar();
subplot(2,3,4);
imagesc(xcat_db.xcat_density(:,:,slice));
title(sprintf('Xcat Vol Density'));
pbaspect([1,1,1]); colorbar();
subplot(2,3,2);
imagesc(final_0_db.final_density(:,:,slice));
title(sprintf('Per Bin Reconstruction Final'));
pbaspect([1,1,1]); colorbar();
subplot(2,3,3);
imagesc(final_1_db.final_density(:,:,slice));
title(sprintf('Accumulation Reconstruction Final'));
pbaspect([1,1,1]); colorbar();
subplot(2,3,5);
imagesc((final_0_db.final_density(:,:,slice)-xcat_db.xcat_density(:,:,slice)));
title(sprintf('Per Bin Reconstruction Error'));
pbaspect([1,1,1]); colorbar();
subplot(2,3,6);
imagesc((final_1_db.final_density(:,:,slice)-xcat_db.xcat_density(:,:,slice)));
title(sprintf('Accumulation Reconstruction Error'));
pbaspect([1,1,1]); colorbar();

%%
slice_to_plot = 5;
relevant_z_indices = [1, 8, 14, 19, 26];
z_names = {'H','O','P','Ca','I'};

figure;
for ii=1:4
    subplot(4,4,ii+0);
    imagesc(squeeze(xcat_db.xcat_fractional_mass(relevant_z_indices(ii),:,:,slice_to_plot)));
    title(['XCAT - ',z_names{ii}]); colorbar();
    subplot(4,4,ii+4);
    imagesc(squeeze(init_db.init_vol_fractional_mass(ii,:,:,slice_to_plot)));
    title(['INIT - ',z_names{ii}]); colorbar();
    subplot(4,4,ii+8);
    imagesc(squeeze(final_0_db.final_fractional_mass(ii,:,:,slice_to_plot)));
    title(['PER BIN - ',z_names{ii}]); colorbar();
    subplot(4,4,ii+12);
    imagesc(squeeze(final_1_db.final_fractional_mass(ii,:,:,slice_to_plot)));
    title(['ACCUM - ',z_names{ii}]); colorbar();
    
end

%%
slice_to_plot = 5;


full_Z = [1:13,15:26,53,82];
relevant_z_indices = [1, 8, 14, 19]; % No Iodine: , 26];
Z = full_Z(relevant_z_indices);
E = 70e-3;%in MeV;
mac = PhotonAttenuationQ(Z, E, 'mac')';
mac_slice = repmat(mac,1,size(xcat_db.xcat_fractional_mass,2),size(xcat_db.xcat_fractional_mass,3));


figure;
subplot(2,2,1);
imagesc(squeeze(sum(mac_slice.*permute(repmat(xcat_db.xcat_density(:,:,slice),1,1,1,4),[4,1,2,3]).*xcat_db.xcat_fractional_mass(relevant_z_indices,:,:,slice_to_plot),1)));
title(sprintf('XCAT Attenuation at E = %d [KeV]',E*1e3));
colormap('bone'); colorbar();
pbaspect([1,1,1]);

subplot(2,2,2);
imagesc(squeeze(sum(mac_slice.*permute(repmat(init_db.init_vol_density(:,:,slice),1,1,1,4),[4,1,2,3]).*init_db.init_vol_fractional_mass(:,:,:,slice_to_plot),1)));
title(sprintf('Linear Inir Attenuation at E = %d [KeV]',E*1e3));
colormap('bone'); colorbar();
pbaspect([1,1,1]);

subplot(2,2,3);
imagesc(squeeze(sum(mac_slice.*permute(repmat(final_0_db.final_density(:,:,slice),1,1,1,4),[4,1,2,3]).*final_0_db.final_fractional_mass(:,:,:,slice_to_plot),1)));
title(sprintf('Final Per-Bin Attenuation at E = %d [KeV]',E*1e3));
colormap('bone'); colorbar();
pbaspect([1,1,1]);

subplot(2,2,4);
imagesc(squeeze(sum(mac_slice.*permute(repmat(final_1_db.final_density(:,:,slice),1,1,1,4),[4,1,2,3]).*final_1_db.final_fractional_mass(:,:,:,slice_to_plot),1)));
title(sprintf('Final Accumulate Attenuation at E = %d [KeV]',E*1e3));
colormap('bone'); colorbar();
pbaspect([1,1,1]);
