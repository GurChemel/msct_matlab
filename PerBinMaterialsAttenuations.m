function [linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names)

addpath('PhotonAttenuation\');
Materials = load('Materials.mat');
load('Materials.mat','Row_names','Element_Z_Count');
Element_names = Row_names(1:(end-1));

active_elements = {'H','O','P','Ca','I'};
num_active_elements = length(active_elements);
relevant_z_indices = zeros(num_active_elements,1);
for ii=1:num_active_elements
    relevant_z_indices(ii) = find(strcmp(Element_names,active_elements{ii}));
end
relevant_z_count = Element_Z_Count(relevant_z_indices);

num_materials = length(Materials_names);
density_vec = zeros(1,num_materials);
fracs_mat   = zeros(num_active_elements,num_materials);
for ii=1:num_materials
    fracs_mat(:,ii) = Materials.(Materials_names{ii})(relevant_z_indices);
    density_vec(ii) = Materials.(Materials_names{ii})(end);
end
density_mat = repmat(density_vec,num_active_elements,1).*fracs_mat;

E_old     = [31.31,41.54,50.39,58.94,69.08,90.85]*1e-3;    %in MeV;
E_src_his = [31.7776   41.6857   50.4728   58.9247   69.1533   90.8185]*1e-3;    %in MeV;
E_dst_his = [33.7091   44.6427   53.5572   61.3088   71.9667   92.8770]*1e-3;    %in MeV;
E_src_exp = [35.5643   48.0523   57.3080   65.7904   78.3154   97.6210]*1e-3;    %in MeV;
E_det_exp = [36.7648   50.4728   59.7230   68.4902   80.8639   98.9679]*1e-3;    %in MeV;
E = E_det_exp;

mac = PhotonAttenuationQ(relevant_z_count, E, 'mac')';
linear_atten_mat = squeeze(sum(repmat(density_mat,1,1,length(E)).*permute(repmat(mac,1,1,num_materials),[1,3,2]),1));

end