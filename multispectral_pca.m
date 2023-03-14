clear; clc;

load('reconstruction_atten_per_bin1.mat')

%%

for slice = 3 %[1,3,6,10]
    per_bin_image(:,1) = reshape(reconstruct_bin_0(:,:,slice),[],1);
    per_bin_image(:,2) = reshape(reconstruct_bin_1(:,:,slice),[],1);
    per_bin_image(:,3) = reshape(reconstruct_bin_2(:,:,slice),[],1);
    per_bin_image(:,4) = reshape(reconstruct_bin_3(:,:,slice),[],1);
    per_bin_image(:,5) = reshape(reconstruct_bin_4(:,:,slice),[],1);
    per_bin_image(:,6) = reshape(reconstruct_bin_5(:,:,slice),[],1);

    figure;
    imagesc(per_bin_image);
    colormap('bone'); pbaspect([1,1,1]);
    
end

%%
[coeffs, a, b] = pca(per_bin_image.');

var_per_bin = var(coeffs);
power_per_bin = sum(abs(coeffs).^2)/3600;

norm_var = var_per_bin./power_per_bin;

%%
figure;
for ii=1:5
    subplot(2,3,ii)
    imagesc(reshape(coeffs(:,ii),60,60));
    colormap('bone'); pbaspect([1,1,1]);
end
%%

gParams;
energy_centers = [31.3090   41.5364   50.3864   58.9392   69.0782   90.8517];

all_bins_options = nchoosek(1:6,2);
num_options = size(all_bins_options,1);

rcond_vec = zeros(num_options,1);
for ii=1:num_options
    selected_bins = all_bins_options(ii,:);


    energy_low  = energy_centers(selected_bins(1));
    energy_high = energy_centers(selected_bins(2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting Attenuation Matrix:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mu_water_e_low  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_low);
    mu_water_e_high = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_high);
    mu_bone_e_low   = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_low);
    mu_bone_e_high  = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_high);

    mu_matrix = [mu_water_e_low,mu_bone_e_low;mu_water_e_high,mu_bone_e_high]
    
    rcond_vec(ii) = rcond(mu_matrix);

end

fprintf('Bin %d and %d. Rcond: %.4f\n',[all_bins_options , rcond_vec].');

% plot(rcond_vec)

%%
gParams;
figure;

energy_to_plot_vec = (30:0.5:120)';

blood_att_vec  = interp1((1e3)*blood_atten(:,1),blood_atten(:,3),energy_to_plot_vec);
water_att_vec  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_to_plot_vec);
muscle_att_vec = interp1((1e3)*muscle_atten(:,1),muscle_atten(:,3),energy_to_plot_vec);
bone_att_vec   = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_to_plot_vec);
iodine_att_vec = interp1((1e3)*iodine_atten(:,1),iodine_atten(:,3),energy_to_plot_vec);

iodine_att_vec=bone_att_vec;

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

gParams;
energy_centers = [31.3090   41.5364   50.3864   58.9392   69.0782   90.8517];

all_bins_options = nchoosek(1:6,3);
num_options = size(all_bins_options,1);

rcond_vec = zeros(num_options,1);
for ii=1:num_options
    selected_bins = all_bins_options(ii,:);


    energy_low  = energy_centers(selected_bins(1));
    energy_mid  = energy_centers(selected_bins(2));
    energy_high = energy_centers(selected_bins(3));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting Attenuation Matrix:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mu_blood_e_low   = interp1((1e3)*iodine_atten(:,1),iodine_atten(:,3),energy_low);
    mu_blood_e_mid   = interp1((1e3)*iodine_atten(:,1),iodine_atten(:,3),energy_mid);
    mu_blood_e_high  = interp1((1e3)*iodine_atten(:,1),iodine_atten(:,3),energy_high);
    mu_muscle_e_low  = interp1((1e3)*muscle_atten(:,1),muscle_atten(:,3),energy_low);
    mu_muscle_e_mid  = interp1((1e3)*muscle_atten(:,1),muscle_atten(:,3),energy_mid);
    mu_muscle_e_high = interp1((1e3)*muscle_atten(:,1),muscle_atten(:,3),energy_high);
    mu_bone_e_low    = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_low);
    mu_bone_e_mid    = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_mid);
    mu_bone_e_high   = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_high);

    mu_matrix = [mu_blood_e_low,mu_muscle_e_low,mu_bone_e_low;mu_blood_e_mid,mu_muscle_e_mid,mu_bone_e_mid;mu_blood_e_high,mu_muscle_e_high,mu_bone_e_high];

    inv_mu_matrix = mu_matrix^-1;
    
    rcond_vec(ii) = rcond(mu_matrix);

end

fprintf('Bins (%d,%d,%d). Rcond: %.7f\n',[all_bins_options , rcond_vec].');
