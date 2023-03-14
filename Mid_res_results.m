clear; clc;

result{1} = load('xcat_reduced.mat');
result{2} = load('Single_0init_vol.mat');
result{3} = load('Single_0final_vol.mat');
result{4} = load('single_1final_vol.mat');

% %%
% DualSliderImshow(result{1}.xcat_density,result{2}.init_vol_density)
%%
init_combined_images = cat(2,result{1}.xcat_density,result{2}.init_vol_density);
% combined_images = cat(2,result{1}.init_vol_density,result{2}.init_vol_density);

num_to_plot = 12;
num_rows = 3;
indices = round(linspace(1,45,num_to_plot));
figure;
for ii=1:num_to_plot
    subplot(num_rows,num_to_plot/num_rows,ii);
    imagesc(init_combined_images(:,:,indices(ii)));
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    title(['Slice ',num2str(indices(ii)),' From 90']);
end

SliderImshow(init_combined_images);

%%
final_combined_images = cat(2,result{1}.xcat_density,result{2}.init_vol_density,result{3}.final_density,result{4}.final_density);
% combined_images = cat(2,result{1}.init_vol_density,result{2}.init_vol_density);

num_to_plot = 1; %3;
num_rows = 1; %3;
indices = 20; % round(linspace(10,30,num_to_plot));
figure;
for ii=1:num_to_plot
    subplot(num_rows,num_to_plot/num_rows,ii);
    imagesc(final_combined_images(:,:,indices(ii)));
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    title(['Slice ',num2str(indices(ii)),' From 90']);
end
%%
indices_to_plot = [1,2,3,4];
% titles_vec = {'H', 'C', 'N', 'O', 'P', 'Ca'};
titles_vec = {'H', 'O', 'P', 'Ca'};
num_to_plot = length(indices_to_plot);
slice_to_plot = 20;
for ii=1:num_to_plot
    subplot(2,num_to_plot,ii);
    imshow(squeeze(result{3}.final_fractional_mass(indices_to_plot(ii),:,:,slice_to_plot)).*result{3}.final_density(:,:,slice_to_plot),[])
    title(titles_vec{indices_to_plot(ii)})
    subplot(2,num_to_plot,ii+num_to_plot);
    imshow(squeeze(result{4}.final_fractional_mass(indices_to_plot(ii),:,:,slice_to_plot)).*result{4}.final_density(:,:,slice_to_plot),[])
end

%%
slice_to_plot = 20;
[dect_water, dect_bone] = Water_Bone_Decomp_fun(slice_to_plot);

multispectral_water = squeeze(result{3}.final_fractional_mass(1,:,:,slice_to_plot)).*result{3}.final_density(:,:,slice_to_plot);
multispectral_bone = squeeze(result{3}.final_fractional_mass(4,:,:,slice_to_plot)).*result{3}.final_density(:,:,slice_to_plot);

bone_max_factor = 0.6;
water_max_factor = 0.14;

multispectral_bone_avg = multispectral_bone / mean(multispectral_bone,'all');
multispectral_water_avg = multispectral_water / mean(multispectral_water,'all');

figure;
if (isempty(water_max_factor))
    subplot(2,3,1); imagesc(dect_water); pbaspect([1,1,1]); title('DECT water'); colorbar();
    subplot(2,3,2); imagesc(multispectral_water); pbaspect([1,1,1]); title('Multispectral water'); colorbar();
else
    subplot(2,3,1); imagesc(dect_water,[0,water_max_factor]); pbaspect([1,1,1]); title('DECT water'); colorbar();
    subplot(2,3,2); imagesc(multispectral_water,[0,water_max_factor]); pbaspect([1,1,1]); title('Multispectral water'); colorbar();
end
if (isempty(bone_max_factor))
    subplot(2,3,4); imagesc(dect_bone); pbaspect([1,1,1]); title('DECT bone'); colorbar();
    subplot(2,3,5); imagesc(multispectral_bone); pbaspect([1,1,1]); title('Multispectral bone'); colorbar();
else
    subplot(2,3,4); imagesc(dect_bone,[0,bone_max_factor]); pbaspect([1,1,1]); title('DECT bone'); colorbar();
    subplot(2,3,5); imagesc(multispectral_bone,[0,bone_max_factor]); pbaspect([1,1,1]); title('Multispectral bone'); colorbar();
end

H2O_ratio = 1;0.126;
P2Ca_ratio = 0.45;

water_gt = squeeze(result{1}.xcat_fractional_mass(1,:,:,slice_to_plot)*H2O_ratio+result{1}.xcat_fractional_mass(8,:,:,slice_to_plot)*(1-H2O_ratio));
bone_gt  = squeeze(result{1}.xcat_fractional_mass(14,:,:,slice_to_plot)*P2Ca_ratio+result{1}.xcat_fractional_mass(19,:,:,slice_to_plot)*(1-P2Ca_ratio));

subplot(2,3,3); imagesc(water_gt); pbaspect([1,1,1]); title('GT water'); colorbar();
subplot(2,3,6); imagesc(bone_gt); pbaspect([1,1,1]); title('GT bone'); colorbar();


%%
% clc;
% 
% slice_number = 45;
% no_scatter_image  = result{1}.init_vol_density(:,:,slice_number);
% yes_scatter_image = result{2}.init_vol_density(:,:,slice_number);
% max_in_scatter = max(yes_scatter_image,[],'all');
% yes_scatter_image_normalized = yes_scatter_image/max_in_scatter;
% 
% gamma_factor = 0.8;
% gamma_thresh = 0.2;
% yes_scatter_image_normalized_corrected = yes_scatter_image_normalized;
% yes_scatter_image_normalized_corrected(yes_scatter_image_normalized_corrected>gamma_thresh) = yes_scatter_image_normalized_corrected(yes_scatter_image_normalized_corrected>gamma_thresh).^(gamma_factor);
% yes_scatter_image_corrected = yes_scatter_image_normalized_corrected*max_in_scatter;
% 
% fprintf('Original Error: %4.2f\n',sum(abs(no_scatter_image-yes_scatter_image),'all'));
% fprintf('Corrected Error: %4.2f\n',sum(abs(no_scatter_image-yes_scatter_image_corrected),'all'));
% 
% combined_image = cat(2,yes_scatter_image,no_scatter_image,yes_scatter_image_corrected);
% imagesc(combined_image)
% title('[Scatter Original] [No Scatter] [Scatter Corrected]')
