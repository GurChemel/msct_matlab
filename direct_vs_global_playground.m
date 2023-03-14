clear;clc;

load('210903/GT_for_Global_VS_Direct.mat');
%%

figure;
slice_number = 1;
x_vec = 1:400;
subplot(6,2,1); imagesc(squeeze(I_Direct(slice_number,x_vec,:)).');
subplot(6,2,3); imagesc(squeeze(I_Direct_and_Global(slice_number,x_vec,:)).');
subplot(6,2,5); imagesc((squeeze(I_Direct_and_Global(slice_number,x_vec,:))-squeeze(I_Direct(slice_number,x_vec,:))).');
slice_number = 4;
x_vec = 501:900;
subplot(6,2,2); imagesc(squeeze(I_Direct(slice_number,x_vec,:)).');
subplot(6,2,4); imagesc(squeeze(I_Direct_and_Global(slice_number,x_vec,:)).');
subplot(6,2,6); imagesc((squeeze(I_Direct_and_Global(slice_number,x_vec,:))-squeeze(I_Direct(slice_number,x_vec,:))).');

% figure;
slice_number = 101;
x_vec = 501:900;
subplot(6,2,6+1); imagesc(squeeze(I_Direct(slice_number,x_vec,:)).');
subplot(6,2,6+3); imagesc(squeeze(I_Direct_and_Global(slice_number,x_vec,:)).');
subplot(6,2,6+5); imagesc((squeeze(I_Direct_and_Global(slice_number,x_vec,:))-squeeze(I_Direct(slice_number,x_vec,:))).');
slice_number = 44;
x_vec = 701:1100;
subplot(6,2,6+2); imagesc(squeeze(I_Direct(slice_number,x_vec,:)).');
subplot(6,2,6+4); imagesc(squeeze(I_Direct_and_Global(slice_number,x_vec,:)).');
subplot(6,2,6+6); imagesc((squeeze(I_Direct_and_Global(slice_number,x_vec,:))-squeeze(I_Direct(slice_number,x_vec,:))).');
%%

slice_number = 1;
x_vec = 30:280;

direct_image = squeeze(I_Direct(slice_number,x_vec,:)).';
direct_and_global_image = squeeze(I_Direct_and_Global(slice_number,x_vec,:)).';
global_image = direct_and_global_image-direct_image;

min_val = min([direct_image,direct_and_global_image,global_image],[],'all');
max_val = max([direct_image,direct_and_global_image,global_image],[],'all');

direct_image_norm = imadjust(direct_image,[min_val, max_val],[1,0]);
direct_and_global_image_norm = imadjust(direct_and_global_image,[min_val, max_val],[1,0]);
global_image_norm = direct_image_norm-direct_and_global_image_norm;

write_im_for_pres(direct_image_norm,'Only Direct Photons',2)
write_im_for_pres(direct_and_global_image_norm,'Direct and Global Photons',2)
write_im_for_pres(global_image_norm,'Only Global Photons',2)

x_vec = 1:500;

direct_image = squeeze(I_Direct(slice_number,x_vec,:)).';
direct_and_global_image = squeeze(I_Direct_and_Global(slice_number,x_vec,:)).';
global_image = direct_and_global_image-direct_image;

min_val = min([direct_image,direct_and_global_image,global_image],[],'all');
max_val = max([direct_image,direct_and_global_image,global_image],[],'all');

direct_image_norm = imadjust(direct_image,[min_val, max_val],[1,0]);
direct_and_global_image_norm = imadjust(direct_and_global_image,[min_val, max_val],[1,0]);
global_image_norm = direct_image_norm-direct_and_global_image_norm;

write_im_for_pres(direct_image_norm,'Bigger Only Direct Photons',2)
write_im_for_pres(direct_and_global_image_norm,'Bigger Direct and Global Photons',2)
write_im_for_pres(global_image_norm,'Bigger Only Global Photons',2)

%%

slice_number = 1;

x_lims_mat = [ 1, 1508;...
              30,  280;...
               1,  650];

% figure
for ii=1:size(x_lims_mat,1)
    x_vec = x_lims_mat(ii,1):x_lims_mat(ii,2);

    direct_image = squeeze(I_Direct(slice_number,x_vec,:)).';
    direct_and_global_image = squeeze(I_Direct_and_Global(slice_number,x_vec,:)).';
    global_image = direct_and_global_image-direct_image;

%     subplot(size(x_lims_mat,1),1,ii)
%     imshow(global_image,[]); colormap('bone');
    
    fprintf('Lim: [%d,%d]. E_global/E_tot = %.2f.\n',x_lims_mat(ii,1),x_lims_mat(ii,2),sum(global_image,'all')/sum(direct_and_global_image,'all'))
%     title(sprintf('E_{global}/E_{tot} = %.2f',sum(global_image,'all')/sum(direct_and_global_image,'all')))
end

%%
slice_number = 1;

figure;
direct_image = squeeze(I_Direct(slice_number,:,:)).';
direct_and_global_image = squeeze(I_Direct_and_Global(slice_number,:,:)).';
global_image = direct_and_global_image-direct_image;

shift_val = 222;
shifted_global_image = [global_image(:,(end-shift_val+1):end),global_image(:,1:(end-shift_val))];

limits_image = zeros(size(global_image));
limits_image(:, 30+shift_val) = max(global_image,[],'all');
limits_image(:,280+shift_val) = max(global_image,[],'all');
limits_image(:,754) = max(global_image,[],'all');

imshow(shifted_global_image+limits_image,[]); colormap('bone');
write_im_for_pres(shifted_global_image+limits_image,'All Global',1)
