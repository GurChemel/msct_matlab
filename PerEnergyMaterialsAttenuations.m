function [linear_atten_mat, density_vec, mac] = PerEnergyMaterialsAttenuations(Materials_names, Energies_kev)

addpath('PhotonAttenuation\');
Materials = load('Materials.mat');
load('Materials.mat','Row_names','Element_Z_Count');
Element_names = Row_names(1:(end-1));

active_elements = {'H','C','O','P','Ca'};
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

E = Energies_kev*1e-3;    %in MeV;

mac = PhotonAttenuationQ(relevant_z_count, E, 'mac')';
linear_atten_mat = squeeze(sum(repmat(density_mat,1,1,length(E)).*permute(repmat(mac,1,1,num_materials),[1,3,2]),1));

end