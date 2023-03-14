clear;clc;
load('221027\init_vol.mat')

slice = 10;
%%
no_asg_init_vol_density_norm = no_asg_init_vol_density/max(no_asg_init_vol_density,[],'all');
figure;
for ii=1:4
    subplot(2,2,ii)
    imshow(squeeze(no_asg_init_vol_fractional_mass(ii,:,:,slice)).*no_asg_init_vol_density_norm(:,:,slice))
end

%%
clear;clc;
load('221108\xcat_reduced.mat')

slice = 1;
%%
xcat_density_norm = xcat_density - min(xcat_density,[],'all');
xcat_density_norm = xcat_density_norm/max(xcat_density_norm,[],'all');

idx = [1,6,8,19];
el_name = {'H','C','O','Ca'};
figure;
for ii=1:4
    subplot(2,2,ii)
    im_to_plot = squeeze(xcat_fractional_mass(idx(ii),:,:,slice)).*xcat_density_norm(:,:,slice);
    imshow(im_to_plot); title(el_name{ii});
    imwrite(im_to_plot,['ImagesForPresentations/GT_density_for_',el_name{ii},'.png']);
end


