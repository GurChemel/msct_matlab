clear; clc;
addpath('C:\Users\gchem\Downloads\GCMex');

% Artificial Noise:
load('clean_att_slices.mat');
labels_1 = GraphCutWrapper(attenuation_slices_noised);

% Real Noise:
load('220201\reconstruction_atten_anti_grid_Enabled.mat');
slice = 5;
attenuation_slices_noised(:,:,1) = reconstruct_bin_0(:,:,slice);
attenuation_slices_noised(:,:,2) = reconstruct_bin_1(:,:,slice);
attenuation_slices_noised(:,:,3) = reconstruct_bin_2(:,:,slice);
attenuation_slices_noised(:,:,4) = reconstruct_bin_3(:,:,slice);
attenuation_slices_noised(:,:,5) = reconstruct_bin_4(:,:,slice);
attenuation_slices_noised(:,:,6) = reconstruct_bin_5(:,:,slice);

load('220201\xcat_reduced.mat');

slice = 5;
xcat_id_orig = xcat_id(:,:,slice);

labels_2 = GraphCutWrapper(attenuation_slices_noised);

index_to_id_dict = [23,    1,      2           38             31         18              6];
xcat_id_recovered_1 = index_to_id_dict(labels_1+1);
xcat_id_recovered_2 = index_to_id_dict(labels_2+1);

figure;
subplot(1,3,1)
imagesc(xcat_id_recovered_1);
title('Recovered 1')
subplot(1,3,2)
imagesc(xcat_id_recovered_2);
title('Recovered 2')
subplot(1,3,3)
imagesc(xcat_id_orig);
title('Original')

%%
clear; clc;
addpath('GCMex\');

load('220426\reconstruction_atten_anti_grid_Disabled.mat')

smothness_factor = 0.03;
labels_3d = GraphCut3dPerBinWrapper(reconstruct_per_bin, smothness_factor);

index_to_id_dict = [23,    1,      2           38             31         18              6];
xcat_id_recovered_3d = index_to_id_dict(labels_3d);
%%
figure;
subplot(2,2,1); imshow(reconstruct_per_bin(:,:,5,4),[]);
subplot(2,2,2); imshow(reconstruct_per_bin(:,:,10,4),[]);
subplot(2,2,3); imagesc(labels_3d(:,:,5));
subplot(2,2,4); imagesc(labels_3d(:,:,10));
