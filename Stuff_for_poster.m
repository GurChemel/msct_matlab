clear; clc;

load('xcat_reduced.mat')
% load('C:\Users\gchem\Geant4\XRayToF_Git\Parsing\High_Res_Recovery\xcat_reduced.mat')

subplot(2,2,1); imagesc(squeeze(sum(xcat_density,3))); colormap('bone');
subplot(2,2,2); imagesc(squeeze(xcat_density(:,:,2))); colormap('bone');
subplot(2,2,3); imagesc(squeeze(xcat_density(:,:,5))); colormap('bone');
subplot(2,2,4); imagesc(squeeze(xcat_density(:,:,8))); colormap('bone');

%%
clear; clc;

load('xcat_reduced.mat')

LoadGParams;

% gParamStruct.dict_mat_to_dens_xcat_id(9) = 18;
% gParamStruct.dict_mat_to_dens_names{end+1} = 'Iodine Blood';

id_vec = unique(xcat_id);
id_vec(id_vec==23) = [];

slice = 5;
colors_mat = [[0.3010 0.7450 0.9330]; ...
[0.8500 0.3250 0.0980]; ...	'#D95319'	
[0.9290 0.6940 0.1250]; ...	'#EDB120'	
[0.4940 0.1840 0.5560]; ...	'#7E2F8E'	
[0.4660 0.6740 0.1880]; ...	'#4DBEEE'	
[0.6350 0.0780 0.1840]];

zoom_factor = 4;

% figure;
image_index_vec = [4,5,6,10,11,12];
all_mat_image = zeros(60,60,3);
for ii=1:length(id_vec)
%     subplot(2,6,image_index_vec(ii));
    color_image = repmat(permute(colors_mat(ii,:),[1,3,2]),60,60,1);
    material_colord_image = color_image.*repmat(xcat_id(:,:,slice)==id_vec(ii),1,1,3);
    all_mat_image = all_mat_image + material_colord_image;
    material_colord_image(material_colord_image==0)=0.5;
%     imshow(material_colord_image);
    zoomed_image = zeros(size(material_colord_image,1)*zoom_factor,size(material_colord_image,2)*zoom_factor,3);
    for rgb_ind=1:3
        zoomed_image(:,:,rgb_ind) = kron(material_colord_image(:,:,rgb_ind),ones(zoom_factor));
    end
    imwrite(zoomed_image,sprintf('ImagesForPresentations\\mat_decomp\\material_%d.png',id_vec(ii)))
%     title(gParamStruct.dict_mat_to_dens_names{(gParamStruct.dict_mat_to_dens_xcat_id==id_vec(ii))})
end

zoom_factor = 10;
all_mat_image(all_mat_image==0)=0.5;
zoomed_image = zeros(size(all_mat_image,1)*zoom_factor,size(all_mat_image,2)*zoom_factor,3);
for rgb_ind=1:3
    zoomed_image(:,:,rgb_ind) = kron(all_mat_image(:,:,rgb_ind),ones(zoom_factor));
end
imwrite(zoomed_image,sprintf('ImagesForPresentations\\mat_decomp\\all_materials.png'))
% subplot(1,2,1);
% imshow(all_mat_image);

%%

for ii=1:10
    image_norm = xcat_density(:,:,ii)/max(xcat_density,[],'all');
    dicomwrite(uint16(image_norm*((2^16)-1)),sprintf('test_%d.dcm',ii))
end