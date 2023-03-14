clear; clc;
load('reconstruction_atten_per_bin0.mat')
load('fake_reconstruction.mat')

slice = 5;

imshow([reconstruct_bin_0(:,:,slice),reconstruct_bin_1(:,:,slice),reconstruct_bin_2(:,:,slice),reconstruct_bin_3(:,:,slice),reconstruct_bin_4(:,:,slice),reconstruct_bin_5(:,:,slice);...
        attenuation_slices(:,:,1),attenuation_slices(:,:,2),attenuation_slices(:,:,3),attenuation_slices(:,:,4),attenuation_slices(:,:,5),attenuation_slices(:,:,6)]);
%%
psnr_vec(1) = psnr(double(reconstruct_bin_0(:,:,slice)),attenuation_slices(:,:,1));
psnr_vec(2) = psnr(double(reconstruct_bin_1(:,:,slice)),attenuation_slices(:,:,2));
psnr_vec(3) = psnr(double(reconstruct_bin_2(:,:,slice)),attenuation_slices(:,:,3));
psnr_vec(4) = psnr(double(reconstruct_bin_3(:,:,slice)),attenuation_slices(:,:,4));
psnr_vec(5) = psnr(double(reconstruct_bin_4(:,:,slice)),attenuation_slices(:,:,5));
psnr_vec(6) = psnr(double(reconstruct_bin_5(:,:,slice)),attenuation_slices(:,:,6));
    