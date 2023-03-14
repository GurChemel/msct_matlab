clear; clc;

rec_80  = load('221204\reconstruction_full_spectrum_80kvp.mat');
rec_140 = load('221204\reconstruction_full_spectrum.mat');

rec_multi = load('221204\reconstruction_atten_per_bin0.mat');

slice = 38;

gParams;
%%
rec_80_normed = rec_80.reconstruct_no_asg;
rec_80_normed = rec_80_normed/max(rec_80_normed,[],'all');
rec_80_normed_gamma = rec_80_normed.^0.5;

rec_140_normed = rec_140.reconstruct_no_asg;
rec_140_normed = rec_140_normed/max(rec_140_normed,[],'all');
rec_140_normed_gamma = rec_140_normed.^0.5;

figure;
subplot(1,2,1);
imshow(squeeze(rec_80_normed_gamma(:,slice,:))')

subplot(1,2,2);
imshow(squeeze(rec_140_normed_gamma(:,slice,:))')

%%

reconstruct_high_freq = squeeze(rec_140.reconstruct_no_asg(:,slice,:))';
reconstruct_low_freq  = squeeze(rec_80.reconstruct_no_asg(:,slice,:))';

energy_low  = 50;
energy_high = 60;

% energy_vec_80 = csvread('spectrum_80kev.txt');
% energy_vec_140 = csvread('spectrum_140kev.txt');
% energy_low  = sum(energy_vec_80'.*(1:150))/sum(energy_vec_80); % In KeV;
% energy_high = sum(energy_vec_140'.*(1:150))/sum(energy_vec_140); % In KeV;

mu_water_e_low  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_low);
mu_water_e_high = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_high);
mu_bone_e_low   = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_low);
mu_bone_e_high  = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_high);

mu_matrix = [mu_water_e_low,mu_bone_e_low;mu_water_e_high,mu_bone_e_high];
% mu_matrix = [mu_water_e_low_calculated,mu_bone_e_low_calculated;mu_water_e_high_calculated,mu_bone_e_high_calculated];

inv_mu_matrix = mu_matrix^-1;

a_mat = inv_mu_matrix*[reconstruct_low_freq(:), reconstruct_high_freq(:)].';
a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);


figure;
subplot(2,2,1); imshow(reconstruct_high_freq,[0,1]); title('High Energy Reconstruction');
subplot(2,2,3); imshow(reconstruct_low_freq,[0,1]); title('Low Energy Reconstruction');
subplot(2,2,2); imshow(norm_a1,[0,1]); title('Water Coefficient');
subplot(2,2,4); imshow(norm_a2,[0,1]); title('Bone Coefficient');


imwrite(kron(reconstruct_high_freq,ones(4)),"ImagesForPresentations/dect_reconstruct_high_freq.jpg");
imwrite(kron(reconstruct_low_freq,ones(4)),"ImagesForPresentations/dect_reconstruct_low_freq.jpg");
imwrite(kron(norm_a1,ones(4)),"ImagesForPresentations/dect_norm_a1.jpg");
imwrite(kron(norm_a2,ones(4)),"ImagesForPresentations/dect_norm_a2.jpg");

%%
reconstruct_high_freq = squeeze(rec_multi.reconstruct_bin_3_4_5(:,slice,:))';
reconstruct_low_freq  = squeeze(rec_multi.reconstruct_bin_0_1_2(:,slice,:))';

energy_low  = 50; % 54.5235
energy_high = 75; % 74.1622

% energy_vec = csvread('spectrum.txt');
% energy_low  = sum(energy_vec(45:64)'.*(45:64))/sum(energy_vec(45:64)); % In KeV;
% energy_high = sum(energy_vec(64:89)'.*(64:89))/sum(energy_vec(64:89)); % In KeV;

mu_water_e_low  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_low);
mu_water_e_high = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_high);
mu_bone_e_low   = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_low);
mu_bone_e_high  = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_high);

mu_matrix = [mu_water_e_low,mu_bone_e_low;mu_water_e_high,mu_bone_e_high];

inv_mu_matrix = mu_matrix^-1;

a_mat = inv_mu_matrix*[reconstruct_low_freq(:), reconstruct_high_freq(:)].';
a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));


norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);


figure;
subplot(2,2,1); imshow(reconstruct_high_freq,[0,1]); title('High Energy Reconstruction');
subplot(2,2,3); imshow(reconstruct_low_freq,[0,1]); title('Low Energy Reconstruction');
subplot(2,2,2); imshow(norm_a1,[0,1]); title('Water Coefficient');
subplot(2,2,4); imshow(norm_a2,[0,1]); title('Bone Coefficient');

imwrite(kron(reconstruct_high_freq,ones(4)),"ImagesForPresentations/multi_reconstruct_high_freq.jpg");
imwrite(kron(reconstruct_low_freq,ones(4)),"ImagesForPresentations/multi_reconstruct_low_freq.jpg");
imwrite(kron(norm_a1,ones(4)),"ImagesForPresentations/multi_norm_a1.jpg");
imwrite(kron(norm_a2,ones(4)),"ImagesForPresentations/multi_norm_a2.jpg");

%%
energy_80_vec = csvread('spectrum_80kev.txt');
energy_140_vec = csvread('spectrum_140kev.txt');
plot([energy_80_vec, energy_140_vec])
xlim([0,140])