clear; clc;

% load('reconstruction_atten_per_bin0.mat')
load('C:\Users\gchem\Geant4\XRayToF_Git\Parsing\Knee_Low_Res\211005\reconstruction_atten_anti_grid_Enabled.mat')

%%
energy_centers = [31.3090   41.5364   50.3864   58.9392   69.0782   90.8517];

for slice = 5 %[1,3,6,10]
    per_bin_image(:,:,1) = reconstruct_bin_0(:,:,slice);
    per_bin_image(:,:,2) = reconstruct_bin_1(:,:,slice);
    per_bin_image(:,:,3) = reconstruct_bin_2(:,:,slice);
    per_bin_image(:,:,4) = reconstruct_bin_3(:,:,slice);
    per_bin_image(:,:,5) = reconstruct_bin_4(:,:,slice);
    per_bin_image(:,:,6) = reconstruct_bin_5(:,:,slice);

    figure;
    for ii=1:6
        subplot(2,4,ii);
%         imshow(squeeze(per_bin_image(:,:,ii)),[0,max(per_bin_image,[],'all')]);
        imagesc(squeeze(per_bin_image(:,:,ii)));
        colormap('bone'); pbaspect([1,1,1]);
        title(sprintf('Energy Center: %.2f',energy_centers(ii)))
    end
    subplot(2,4,7);
    imagesc(sum(per_bin_image,3));
    colormap('bone'); pbaspect([1,1,1]);
    title(sprintf('Energy Accumulated'))
    
end


%%

selected_bins = [4,6];

% energy_centers = [26.5, 41.5, 50.5, 59, 70, 98.5];
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
figure;
subplot(2,2,1); 
imagesc(reconstruct_low_freq); title(sprintf('Low Energy - %.2f [KeV]',energy_low));
colormap('bone'); pbaspect([1,1,1]);
subplot(2,2,2); 
imagesc(reconstruct_high_freq); title(sprintf('High Energy - %.2f [KeV]', energy_high));
colormap('bone'); pbaspect([1,1,1]);

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
% subplot(2,2,3);
% imagesc(norm_a1,[0,1]); title('Water Coefficient');
% colormap('bone'); pbaspect([1,1,1]);
% subplot(2,2,4);
% imagesc(norm_a2,[0,1]); title('Bone Coefficient');
% colormap('bone'); pbaspect([1,1,1]);

% bone_idx = medfilt2(norm_a1<norm_a2,[3 3]);
bone_idx = (0.8*norm_a1)<norm_a2;
bone_idx = 0.2<norm_a2;
% bone_idx = medfilt2(bone_idx, [2 2]);
% 
% SE = strel('rectangle',[1 2]);
% bone_idx = imerode(bone_idx,SE);
% bone_idx = imdilate(bone_idx,SE);

subplot(2,2,3);
imagesc(reconstruct_low_freq.*(bone_idx==0)); title('Water Coefficient');
colormap('bone'); pbaspect([1,1,1]);
subplot(2,2,4);
imagesc(reconstruct_low_freq.*(bone_idx==1)); title('Bone Coefficient');
colormap('bone'); pbaspect([1,1,1]);

%%

selected_bins = [1,2,3,4,5,6];

energy_centers = [31.3090   41.5364   50.3864   58.9392   69.0782   90.8517].';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Reconstructions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for ii=selected_bins
    subplot(2,length(energy_centers),ii); 
    imagesc(squeeze(per_bin_image(:,:,ii))); title(sprintf('Energy - %.2f [KeV]',energy_centers(ii)));
    colormap('bone'); pbaspect([1,1,1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Attenuation Matrix:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gParams;

mu_water_e  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_centers(selected_bins));
mu_bone_e   = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_centers(selected_bins));

mu_matrix = [mu_water_e,mu_bone_e];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating wieghts:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

energy_vec = csvread('spectrum.txt');
full_water_atten = interp1((1e3)*water_atten(:,1),water_atten(:,3),1:150)';
full_bone_atten = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),1:150)';

bone_dec_energy_vec = energy_vec.*exp(-full_bone_atten);
water_dec_energy_vec = energy_vec.*exp(-full_water_atten);
edges_vec = [16, 37, 46, 55, 63, 77, 120];

for ii=1:6
    input_sum(ii) = sum(energy_vec(edges_vec(ii):edges_vec(ii+1)));
    bone_output_sum(ii) = sum(bone_dec_energy_vec(edges_vec(ii):edges_vec(ii+1)));
    water_output_sum(ii) = sum(water_dec_energy_vec(edges_vec(ii):edges_vec(ii+1)));
end

% atten_data_mat_for_plot = [energy_vec , bone_dec_energy_vec, water_dec_energy_vec];
% figure;
% subplot(1,2,1); plot(atten_data_mat_for_plot);
% title('Spectrum'); legend('Source','Bone Detector','Water Detector');
% subplot(1,2,2); plot(exp(-full_bone_atten)); title('Bone Attenuation');

% weight_vec = exp(-(1.8*mu_bone_e)); %[low_weight;low_weight;low_weight;high_weight;low_weight;high_weight];
weight_vec = (bone_output_sum./input_sum)';
weight_vec = weight_vec(selected_bins);

% mu_matrix = mu_matrix.*(repmat(weight_vec,1,2).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pseudo-Inverse:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inv_mu_matrix = pinv(mu_matrix);
% covariance_matrix = diag(var(reshape(per_bin_image,[],6)));
% mu_matrix(1,:) = [0, 0];
% mu_matrix(2,:) = [0, 0];
% mu_matrix(3,:) = [0, 0];
% mu_matrix(4,:) = [0, 0];
% covariance_matrix = eye(6);
covariance_matrix = diag(1./[761828, 798509, 818354, 895294, 917473, 1009081]);
% inv_mu_matrix = (mu_matrix.'*(covariance_matrix)^-1*mu_matrix)\(mu_matrix.')*((covariance_matrix)^-1);
inv_mu_matrix = (mu_matrix.'*covariance_matrix*mu_matrix)\(mu_matrix.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weight_vec = [0,0,0.1,.8,1,1]';
% a_mat = lscov(mu_matrix,reshape(per_bin_image(:,:,selected_bins),[],length(selected_bins)).',weight_vec);
a_mat = inv_mu_matrix*reshape(per_bin_image(:,:,selected_bins),[],length(selected_bins)).';
a_1 = reshape(a_mat(1,:),size(per_bin_image(:,:,1)));
a_2 = reshape(a_mat(2,:),size(per_bin_image(:,:,1)));

norm_a1 = max(a_1./max([0,max(a_1,[],'all')]),0);
norm_a2 = max(a_2./max([0,max(a_2,[],'all')]),0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Decompositions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,6,9);
imagesc(norm_a1,[0,1]); title('Water Coefficient');
colormap('bone'); pbaspect([1,1,1]);
subplot(2,6,10);
imagesc(norm_a2,[0,1]); title('Bone Coefficient');
colormap('bone'); pbaspect([1,1,1]);