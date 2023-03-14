clear;clc;
% load('xcat_full.mat');
load('210903/xcat_reduced.mat');
active_elements_struct=load('xcat_elements.mat');
active_elements_cell = fieldnames(active_elements_struct);
active_materials_struct=load('xcat_materials.mat');
active_materials_cell = fieldnames(active_materials_struct);

% load('Single_0init_vol.mat');

%%
slices_to_plot = [2,4,6,8];

plot_1_or_save_0 = 0;

if plot_1_or_save_0
    figure;
end

xcat_density_norm = xcat_density/max(xcat_density,[],'all');
xcat_density_norm(xcat_density_norm<0)=0;

for ii=1:length(slices_to_plot)

    slice_density = xcat_density_norm(:,:,slices_to_plot(ii));
    title_str = sprintf('GT Density Slice %d from %d',slices_to_plot(ii),size(xcat_density_norm,3));
    if plot_1_or_save_0
        subplot(2,2,ii); imagesc(slice_density);
        colormap('bone'); colorbar(); pbaspect([1,1,1]);
        title(title_str);
    end
    if plot_1_or_save_0==0
        file_name_str = title_str;
        file_name_str(file_name_str==' ')='_';
        file_name_str = ['ImagesForPresentations/',file_name_str,'.png'];
        imwrite(kron(slice_density,ones(4)),file_name_str);
    end
end


%%

slice_to_plot = 5;
slice_fractional_mass = xcat_fractional_mass(:,:,:,slice_to_plot);
slice_density = xcat_density(:,:,slice_to_plot);


Z = [1:13,15:26,53,82];
E = [50,70,90]*1e-3;%in MeV;
%meac = PhotonAttenuationQ(Z, E, 'meac');
figure;
subplot(2,2,1); imagesc(slice_density); colormap('bone'); colorbar(); title('XCAT Density Slice'); pbaspect([1,1,1]);

for ii=1:3
    mac = PhotonAttenuationQ(Z, E(ii), 'mac')';

    mac_slice = repmat(mac,1,size(slice_density,1),size(slice_density,2));

    slice_per_z_attenuation = mac_slice.*slice_fractional_mass;
    slice_attenuation = squeeze(sum(slice_per_z_attenuation,1));

    subplot(2,2,1+ii); imagesc(slice_attenuation);
    colormap('bone'); colorbar(); pbaspect([1,1,1]);
    title(sprintf('XCAT Linear Attenuation Slice. E = %d [KeV]',E(ii)*1e3));
end

%%
slice_to_plot = 1;
slice_fractional_mass = xcat_fractional_mass(:,:,:,slice_to_plot);
slice_density = xcat_density(:,:,slice_to_plot);

figure
for ii=1:27
   subplot(3,9,ii);
   element_density = squeeze(slice_fractional_mass(ii,:,:)).*slice_density;
   imagesc(element_density);
   colorbar();
   title(active_elements_cell{ii});
end
%%
figure
slice = 1;
for slice_to_plot=[1,3,5,10]
    slice_fractional_mass = xcat_fractional_mass(:,:,:,slice_to_plot);
    slice_density = xcat_density(:,:,slice_to_plot);
    element = 1;
    for element_to_plot=[1,8,14,19]
       subplot(4,4,(slice-1)*4+element);
       element_density = squeeze(slice_fractional_mass(element_to_plot,:,:)).*slice_density;
       imagesc(element_density);
       colorbar();
       title(active_elements_cell{element_to_plot});
       element = element + 1;
    end
    slice = slice + 1;
end

%%
elements_mass = xcat_fractional_mass;%.*repmat(permute(xcat_density,[4,1,2,3]),27,1,1,1);
figure
slice = 1;
for slice_to_plot=[1,3,5,10]
    slice_elements_mass = elements_mass(:,:,:,slice_to_plot);
    element = 1;
    for element_to_plot=[1,8,14,19]
       subplot(4,4,(slice-1)*4+element);
       element_mass = squeeze(slice_elements_mass(element_to_plot,:,:));
       imagesc(element_mass,[0,1]);
       colorbar();
       title(active_elements_cell{element_to_plot});
       element = element + 1;
    end
    slice = slice + 1;
end


%%
uniqe_id_value = unique(xcat_id(:,:,slice_to_plot));
num_uniques = length(uniqe_id_value);
unique_locations = zeros(num_uniques,1);
active_elements_categorical = categorical(active_elements_cell);

figure
for ii=1:num_uniques
    unique_locations(ii) = find(xcat_id(:,:,slice_to_plot)==uniqe_id_value(ii),1);
    subplot(2,4,ii)
    bar(active_elements_categorical,slice_fractional_mass(:,unique_locations(ii)))
    material_name = active_materials_cell{uniqe_id_value(ii)};
    material_name(material_name=='_')=' ';
    title(material_name);
end
%%

slice_id = xcat_id(:,:,slice_to_plot);
uniqe_id_value = unique(slice_id);
num_uniques = length(uniqe_id_value);
simplified_id_map = slice_id;
for ii=1:num_uniques
    simplified_id_map(slice_id==uniqe_id_value(ii))=ii;
    material_name = active_materials_cell{uniqe_id_value(ii)};
    material_name(material_name=='_')=' ';
    simplified_legend{ii} = material_name;
end

figure;
cmap = hsv(num_uniques);
image(simplified_id_map);
colormap(cmap);
hold on;
for K = 1 : num_uniques; hidden_h(K) = surf(uint8(K-[1 1;1 1]), 'edgecolor', 'none'); end
hold off
uistack(hidden_h, 'bottom');
legend(hidden_h, simplified_legend )

%%

xcat_slice = xcat_density(:,:,5);
init_slice = init_vol_density(:,:,5);

figure;
subplot(2,2,1); imagesc(xcat_slice); colorbar(); title('XCAT Slice'); pbaspect([1,1,1]);
subplot(2,2,2); imagesc(init_slice); colorbar(); title('Init Slice'); pbaspect([1,1,1]);
subplot(2,2,3); imagesc((init_slice-xcat_slice)); colorbar(); title('Error Slice'); pbaspect([1,1,1]);
subplot(2,2,4); imagesc(medfilt2((init_slice-xcat_slice))); colorbar(); title('Error Slice Filtered'); pbaspect([1,1,1]);

xcat_density_vec = xcat_slice(:);%+0.01*rand(numel(xcat_slice),1);
init_density_vec = init_slice(:);

max_idx = max([max(xcat_density_vec),max(init_density_vec)]);
%%
% subplot(2,2,4);
figure;
scatter(xcat_density_vec,init_density_vec); hold on;
xlabel('Ground Truth')
ylabel('Linear Initialization')
plot([0,max_idx],[0,max_idx]); hold off;
ylim([0,max_idx]);
xlim([0,max_idx]);
