clear;clc;

load('reconstruction_atten_per_bin0.mat')

slice_to_plot = 5;
per_bin_image(:,:,1) = reconstruct_bin_0(:,:,slice_to_plot);
per_bin_image(:,:,2) = reconstruct_bin_1(:,:,slice_to_plot);
per_bin_image(:,:,3) = reconstruct_bin_2(:,:,slice_to_plot);
per_bin_image(:,:,4) = reconstruct_bin_3(:,:,slice_to_plot);
per_bin_image(:,:,5) = reconstruct_bin_4(:,:,slice_to_plot);
per_bin_image(:,:,6) = reconstruct_bin_5(:,:,slice_to_plot);

figure
for ii=1:6
    subplot(2,3,ii);
    imagesc(per_bin_image(:,:,ii)); colormap('bone');
end

%%
clear;clc;
%%
energy_low  = 50;
energy_high = 80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Reconstructions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slice_to_plot = 5;
single_init = load(['Single_',num2str(energy_low),'init_vol.mat'],'reconstructed');
reconstruct_low_freq = single_init.reconstructed(:,:,slice_to_plot);
single_init = load(['Single_',num2str(energy_high),'init_vol.mat'],'reconstructed');
reconstruct_high_freq = single_init.reconstructed(:,:,slice_to_plot);
clear single_init;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Reconstructions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,4,1); 
title_str = 'MONO Low Energy - 50 [KeV]';
imagesc(reconstruct_low_freq); title(title_str);
colormap('bone'); pbaspect([1,1,1]);
write_im_for_pres(reconstruct_low_freq,title_str)
subplot(2,4,2); 
title_str = 'MONO Low Energy - 80 [KeV]';
imagesc(reconstruct_high_freq); title(title_str);
colormap('bone'); pbaspect([1,1,1]);
write_im_for_pres(reconstruct_high_freq,title_str)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Attenuation Matrix:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gParams;

mu_water_e_low  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_low);
mu_water_e_high = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_high);
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
subplot(2,4,5);
title_str = 'MONO Water Coefficient';
imagesc(norm_a1,[0,1]); title(title_str);
colormap('bone'); pbaspect([1,1,1]);
write_im_for_pres(norm_a1,title_str)

subplot(2,4,6);
title_str = 'MONO Bone Coefficient';
imagesc(norm_a2,[0,1]); title(title_str);
colormap('bone'); pbaspect([1,1,1]);
write_im_for_pres(norm_a2,title_str)

mono_norm_a1 = norm_a1;
mono_norm_a2 = norm_a2;
%%

load('reconstruction_atten_per_bin0.mat')

per_bin_image(:,:,1) = reconstruct_bin_0(:,:,slice_to_plot);
per_bin_image(:,:,2) = reconstruct_bin_1(:,:,slice_to_plot);
per_bin_image(:,:,3) = reconstruct_bin_2(:,:,slice_to_plot);
per_bin_image(:,:,4) = reconstruct_bin_3(:,:,slice_to_plot);
per_bin_image(:,:,5) = reconstruct_bin_4(:,:,slice_to_plot);
per_bin_image(:,:,6) = reconstruct_bin_5(:,:,slice_to_plot);

selected_bins = [4,6];

energy_centers = [31.3090   41.5364   50.3864   58.9392   69.0782   90.8517];

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
subplot(2,4,3); 
title_str = sprintf('MULTI Low Energy - %.2f [KeV]',energy_low);
imagesc(reconstruct_low_freq); title(title_str);
colormap('bone'); pbaspect([1,1,1]);
write_im_for_pres(reconstruct_low_freq,title_str)
subplot(2,4,4); 
title_str = sprintf('MULTI Low Energy - %.2f [KeV]',energy_high);
imagesc(reconstruct_high_freq); title(title_str);
colormap('bone'); pbaspect([1,1,1]);
write_im_for_pres(reconstruct_high_freq,title_str)

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

mu_water = interp1((1e3)*blood_atten(:,1),blood_atten(:,3),energy_centers)';
mu_bone = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_centers)';
weight_vector = [0,0.1,0.5,0.3,0.4,0.5];
% weight_vector = (1./mu_water).^0.7;
% weight_vector = (1./std(reshape(per_bin_image,[],6)));
mu_matrix_full = diag(weight_vector)*[mu_water,mu_bone];
inv_mu_matrix_full = ((mu_matrix_full'*mu_matrix_full)^-1)*mu_matrix_full';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a_mat = inv_mu_matrix*[reconstruct_low_freq(:), reconstruct_high_freq(:)].';
a_mat = inv_mu_matrix_full*(reshape(per_bin_image,[],6)');
a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,4,7);
imagesc(norm_a1,[0,1]); t=title('MULTI Water Coefficient');
colormap('bone'); pbaspect([1,1,1]);
write_im_for_pres(norm_a1,t.String)
subplot(2,4,8);
imagesc(norm_a2,[0,1]); t=title('MULTI Bone Coefficient');
colormap('bone'); pbaspect([1,1,1]);
write_im_for_pres(norm_a2,t.String)

multi_norm_a1 = norm_a1;
multi_norm_a2 = norm_a2;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting MONO Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;

mono_bone_thresh = 0.4;
mono_air_thresh = 0.01;
mono_bone = mono_norm_a1.*(mono_norm_a1>mono_bone_thresh);
mono_water = mono_norm_a1.*(mono_norm_a1>mono_air_thresh).*(mono_norm_a1<=mono_bone_thresh);

subplot(2,2,1);
title_str = 'MONO Water Coefficient';
imagesc(mono_water,[0,1]); title(title_str);
colormap('bone'); pbaspect([1,1,1]);
subplot(2,2,2);
title_str = 'MONO Bone Coefficient';
imagesc(mono_bone,[0,1]); title(title_str);
colormap('bone'); pbaspect([1,1,1]);
%%
multi_bone_thresh = 0.40;
multi_air_thresh = 0.01;
multi_bone = multi_norm_a1.*(multi_norm_a1>multi_bone_thresh);
multi_water = multi_norm_a1.*(multi_norm_a1>multi_air_thresh).*(multi_norm_a1<=multi_bone_thresh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3);
imagesc(multi_water,[0,1]); t=title('MULTI Water Coefficient');
colormap('bone'); pbaspect([1,1,1]);
subplot(2,2,4);
imagesc(multi_bone,[0,1]); t=title('MULTI Bone Coefficient');
colormap('bone'); pbaspect([1,1,1]);

%%
gParams;
figure;

energy_to_plot_vec = (30:0.5:120)';

blood_att_vec  = interp1((1e3)*blood_atten(:,1),blood_atten(:,3),energy_to_plot_vec);
water_att_vec  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_to_plot_vec);
muscle_att_vec = interp1((1e3)*muscle_atten(:,1),muscle_atten(:,3),energy_to_plot_vec);
bone_att_vec   = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_to_plot_vec);
iodine_att_vec = interp1((1e3)*iodine_atten(:,1),iodine_atten(:,3),energy_to_plot_vec);

subplot(2,1,1);
plot(repmat(energy_to_plot_vec,1,5),[blood_att_vec,water_att_vec,muscle_att_vec,bone_att_vec,iodine_att_vec])
legend({'Blood','Water','Muscle','Bone','Iodine'});
xlabel('Energy [KeV]');
ylabel('Linear Attenuation Coeff [cm^2/g]')

subplot(2,1,2);
plot(repmat(energy_to_plot_vec,1,3),[blood_att_vec,water_att_vec,muscle_att_vec])
legend({'Blood','Water','Muscle'});
xlabel('Energy [KeV]');
ylabel('Linear Attenuation Coeff [cm^2/g]')

%%
load('xcat_reduced.mat')

slice_to_plot = 5;
xcat_id_slice = xcat_id(:,:,slice_to_plot);
bone_gt_id = (xcat_id_slice==38);
water_gt_id = ((xcat_id_slice~=23).*(xcat_id_slice~=38))==1;


figure;
subplot(2,3,1); imshow(bone_gt_id);
subplot(2,3,2); imshow(a_2,[]);
subplot(2,3,3); plot([0,1],[0,1]); hold on; plot(a_2(:),bone_gt_id(:),'.');

subplot(2,3,4); imshow(water_gt_id);
subplot(2,3,5); imshow(a_1,[]);
subplot(2,3,6); plot([0,1],[0,1]); hold on; plot(a_1(:),water_gt_id(:),'.');
