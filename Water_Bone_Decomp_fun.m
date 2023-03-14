function [water, bone] = Water_Bone_Decomp_fun(slice_to_plot)
load('reconstruction_atten_low_high.mat','asg_reconstruct_low_freq','asg_reconstruct_high_freq');

reconstruct_high_freq = double(asg_reconstruct_high_freq(:,:,slice_to_plot));
reconstruct_low_freq  = double(asg_reconstruct_low_freq(:,:,slice_to_plot));

gParams;
%%

energy_low  = 53;
energy_high = 65;

% water_mass_density_in_room_temp = 0.99802;
% water_atten(:,2) = water_mass_density_in_room_temp*water_atten(:,2);
% bone_mass_density_in_room_temp = 1.8;
% bone_atten(:,2)  = bone_mass_density_in_room_temp*bone_atten(:,2); % bone_Attem = mu/ro

mu_water_e_low  = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_low);
mu_water_e_high = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_high);
mu_bone_e_low   = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_low);
mu_bone_e_high  = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_high);

mu_matrix = [mu_water_e_low,mu_bone_e_low;mu_water_e_high,mu_bone_e_high];

inv_mu_matrix = mu_matrix^-1;

%%

a_mat = inv_mu_matrix*[reconstruct_low_freq(:), reconstruct_high_freq(:)].';
a_1 = reshape(a_mat(1,:),size(reconstruct_high_freq));
a_2 = reshape(a_mat(2,:),size(reconstruct_high_freq));

%%
water = max(a_1./max([0,max(a_1,[],'all')]),0);
bone = max(a_2./max([0,max(a_2,[],'all')]),0);

end