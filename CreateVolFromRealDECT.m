clear; clc;

source_50 = 'C:\Users\gchem\Downloads\DECT_case_from_rajiv\For_Pydicom\100_kvp';
source_75 = 'C:\Users\gchem\Downloads\DECT_case_from_rajiv\For_Pydicom\150_kvp';
source_50_files = dir([source_50,'\*.dcm']);
info_50 = dicominfo([source_50,'\',source_50_files(1).name]);
source_75_files = dir([source_75,'\*.dcm']);
info_75 = dicominfo([source_75,'\',source_75_files(1).name]);

high_50_kev_hu_images_orig = double(permute(dicomreadVolume(source_50),[1,2,4,3]))*info_50.RescaleSlope+info_50.RescaleIntercept;
high_75_kev_hu_images_orig = double(permute(dicomreadVolume(source_75),[1,2,4,3]))*info_75.RescaleSlope+info_75.RescaleIntercept;

%%

high_50_kev_hu_images = high_50_kev_hu_images_orig;%(1:2:end,1:2:end,:);
high_75_kev_hu_images = high_75_kev_hu_images_orig;%(1:2:end,1:2:end,:);

num_pixels_xy = size(high_50_kev_hu_images,1);

energies = [50,75]; % According to rajiv, E = kVp/2;

high_50_kev_hu_images(high_50_kev_hu_images<-1000)=-1000;
high_75_kev_hu_images(high_75_kev_hu_images<-1000)=-1000;

[mu_water, ~] = PerEnergyMaterialsAttenuations({'Water'}, energies);
[mu_air, ~] = PerEnergyMaterialsAttenuations({'Air'}, energies);

high_50_kev_mu_images = mu_water(1)+((mu_water(1)-mu_air(1))*double(high_50_kev_hu_images)/1000);
high_75_kev_mu_images = mu_water(2)+((mu_water(2)-mu_air(2))*double(high_75_kev_hu_images)/1000);


ct_to_dens = ([[-5000,	0.0];   ...
               [-1000,	0.01];  ...
               [-400,	0.602]; ...
               [-150,	0.924]; ...
               [100,	1.075]; ...
               [300,	1.145]; ...
               [2000,	1.683]; ...
               [4927,   3.379]; ...
               [66000,  7.8]]);
               

%%

volume_slices = 377:400;
num_slices = length(volume_slices);
relevant_high_50_kev_mu_images = high_50_kev_mu_images(:,:,volume_slices);
relevant_high_75_kev_mu_images = high_75_kev_mu_images(:,:,volume_slices);

[mu_mat, ~] = PerEnergyMaterialsAttenuations({'Air','Water','Dryrib'}, energies);

[~,relevant_50_id] = min(abs(repmat(permute(mu_mat(:,1),[4,3,2,1]),num_pixels_xy,num_pixels_xy,num_slices,1)-relevant_high_50_kev_mu_images),[],4);
[~,relevant_75_id] = min(abs(repmat(permute(mu_mat(:,2),[4,3,2,1]),num_pixels_xy,num_pixels_xy,num_slices,1)-relevant_high_75_kev_mu_images),[],4);

relevant_high_50_kev_mu_images(relevant_50_id==1)=0;
relevant_high_75_kev_mu_images(relevant_75_id==1)=0;

[X,Y] = meshgrid(-(num_pixels_xy/2):((num_pixels_xy/2)-1),-(num_pixels_xy/2):((num_pixels_xy/2)-1));

mask = ones(num_pixels_xy,num_pixels_xy,num_slices);
mask(repmat(((((X.^2)+(Y.^2))<125^2)+(Y<0))>0,1,1,num_slices)) = 0;

relevant_high_50_kev_mu_images(mask==1)=0;
relevant_high_75_kev_mu_images(mask==1)=0;

relevant_high_50_kev_mu_images(relevant_high_50_kev_mu_images==0)=Inf;

C_val = relevant_high_75_kev_mu_images./relevant_high_50_kev_mu_images;

num = mu_mat(3,2)-mu_mat(3,1);
denum = C_val*(mu_mat(2,1)-mu_mat(3,1))-(mu_mat(2,2)-mu_mat(3,2));

m_1 = num/denum;

% inv_mu_matrix = (mu_mat(2:3,:)')^-1;
% 
% a_mat = inv_mu_matrix*[relevant_high_50_kev_mu_images(:), relevant_high_75_kev_mu_images(:)].';
% a_1 = reshape(a_mat(1,:),size(relevant_high_75_kev_mu_images));
% a_2 = reshape(a_mat(2,:),size(relevant_high_75_kev_mu_images));
% 
% SliderImshow(cat(2,a_1,a_2))


%%

slice = 10;
relevant_high_50_kev_mu_images = high_50_kev_mu_images(:,:,volume_slices);
relevant_high_75_kev_mu_images = high_75_kev_mu_images(:,:,volume_slices);

[mu_mat, ~] = PerEnergyMaterialsAttenuations({'Air','Water','Dryrib'}, energies);

[X,Y] = meshgrid(-(num_pixels_xy/2):((num_pixels_xy/2)-1),-(num_pixels_xy/2):((num_pixels_xy/2)-1));

mask = zeros(num_pixels_xy,num_pixels_xy);
mask(((((X.^2)+(Y.^2))<125^2)+(Y<0))>0) = 1;

slices_combined = [mask.*relevant_high_50_kev_mu_images(:,:,slice),mask.*relevant_high_75_kev_mu_images(:,:,slice)];

ids_combined = [relevant_50_id(:,:,slice),relevant_75_id(:,:,slice)];

% air_thresh = 0.1;
% water_thresh = 0.18;
% slices_air = slices_combined<air_thresh;
% slices_water = (slices_combined>=air_thresh).*(slices_combined<water_thresh);
% slices_bone = 1-slices_water-slices_air;

figure;
subplot(2,2,1);
imagesc(slices_combined); colormap('bone');
subplot(2,2,2);
imagesc(slices_combined.*(ids_combined==1)); colormap('bone');
subplot(2,2,3);
imagesc(slices_combined.*(ids_combined==2)); colormap('bone');
subplot(2,2,4);
imagesc(slices_combined.*(ids_combined==3)); colormap('bone');

%%

[X,Y] = meshgrid(-(num_pixels_xy/2):((num_pixels_xy/2)-1),-(num_pixels_xy/2):((num_pixels_xy/2)-1));
mask = zeros(num_pixels_xy,num_pixels_xy);
mask(((((X.^2)+(Y.^2))<125^2)+(Y<0))>0) = 1;
mask = repmat(mask,1,1,length(volume_slices));

mask(relevant_50_id==1)=0;
mask(relevant_75_id==1)=0;

slice = 10;
relevant_high_50_kev_mu_images = mask.*high_50_kev_mu_images(:,:,volume_slices);
relevant_high_75_kev_mu_images = mask.*high_75_kev_mu_images(:,:,volume_slices);

[mu_mat, ~] = PerEnergyMaterialsAttenuations({'Air','Water','Dryrib'}, energies);

inv_mu_matrix = (mu_mat(2:3,:)')^-1;

a_mat = inv_mu_matrix*[relevant_high_50_kev_mu_images(:), relevant_high_75_kev_mu_images(:)].';
a_1 = reshape(a_mat(1,:),size(relevant_high_50_kev_mu_images));
a_2 = reshape(a_mat(2,:),size(relevant_high_50_kev_mu_images));

a_1(a_1<0)=0;
a_2(a_2<0)=0;

figure;
subplot(2,4,1); 
imagesc(a_1(:,:,slice)); colormap('bone'); title('water [unitless]'); colorbar;
axis off; pbaspect([1,1,1]);
subplot(2,4,5);
imagesc(a_2(:,:,slice)); colormap('bone'); title('bone [unitless]'); colorbar;
axis off; pbaspect([1,1,1]);

subplot(2,4,2);
imagesc(relevant_high_50_kev_mu_images(:,:,slice)); colormap('bone'); title('low energy orig [1/cm]'); colorbar;
axis off; pbaspect([1,1,1]);
subplot(2,4,6);
imagesc(relevant_high_75_kev_mu_images(:,:,slice)); colormap('bone'); title('high energy orig [1/cm]'); colorbar;
axis off; pbaspect([1,1,1]);

subplot(2,4,3);
imagesc(a_1(:,:,slice)*mu_mat(2,1)+a_2(:,:,slice)*mu_mat(3,1)); colormap('bone'); title('low energy res [1/cm]'); colorbar;
axis off; pbaspect([1,1,1]);
subplot(2,4,7);
imagesc(a_1(:,:,slice)*mu_mat(2,2)+a_2(:,:,slice)*mu_mat(3,2)); colormap('bone'); title('high energy res [1/cm]'); colorbar;
axis off; pbaspect([1,1,1]);

subplot(2,4,4);
imagesc(relevant_high_50_kev_mu_images(:,:,slice)-(a_1(:,:,slice)*mu_mat(2,1)+a_2(:,:,slice)*mu_mat(3,1))); colormap('bone'); title('low energy err [1/cm]'); colorbar;
axis off; pbaspect([1,1,1]);
subplot(2,4,8);
imagesc(relevant_high_75_kev_mu_images(:,:,slice)-(a_1(:,:,slice)*mu_mat(2,2)+a_2(:,:,slice)*mu_mat(3,2))); colormap('bone'); title('high energy err [1/cm]'); colorbar;
axis off; pbaspect([1,1,1]);

%%

water_density = 1;
bone_sensity = 1.55;

full_density = water_density*a_1+bone_sensity*a_2;
water_density = zeros(size(a_1));
bone_density  = zeros(size(a_1));

id_mat = zeros(size(a_1));
id_mat(a_1 > 0.1) = 1;
id_mat(a_2 > 0.1) = 2;
water_density(id_mat==1)=full_density(id_mat==1);
bone_density(id_mat==2) =full_density(id_mat==2);

base_folder = 'dect_phantoms';
dir_path = [base_folder,'/',datestr(date,'YYmmDD'),'/'];
cnt=1;
while exist(dir_path,'dir')
dir_path = [base_folder,'/',datestr(date,'YYmmDD'),'_',num2str(cnt),'/'];
    cnt = cnt + 1;
end
mkdir(dir_path)

for ii=1:size(water_density,3)
    fid = fopen([dir_path,'dect_vol_slice_',sprintf('%03d',ii),'.txt'],'w');
    fprintf(fid,[repmat('%.3f\t',1,num_pixels_xy-1),'%6.3f\n'],water_density(:,:,ii)');
    fclose(fid);
end
for ii=1:size(bone_density,3)
    fid = fopen([dir_path,'dect_vol_slice_',sprintf('%03d',ii+size(water_density,3)),'.txt'],'w');
    fprintf(fid,[repmat('%.3f\t',1,num_pixels_xy-1),'%6.3f\n'],bone_density(:,:,ii)');
    fclose(fid);
end

%%

a_1_sliced = a_1(:,:,slice);
a_2_sliced = a_2(:,:,slice);

water_density_sliced = water_density(:,:,slice);
bone_density_sliced = bone_density(:,:,slice);

fractional = [0.112	0.000
              0.000	0.155
              0.888	0.435
              0.000	0.103
              0.000	0.225];

result = repmat(a_1_sliced,1,1,5).*repmat(permute(fractional(:,1),[2,3,1]),256,256,1)+repmat(a_2_sliced,1,1,5).*repmat(permute(fractional(:,2),[2,3,1]),256,256,1);

result_flat = reshape(result,[],5);

[result_ids, clusters] = kmeans(result_flat,50);

result_quant = reshape(clusters(result_ids,:),256,256,5);

[mu_mat, ~, mac_per_z] = PerEnergyMaterialsAttenuations({'Water','Dryrib'}, 50);

result_attn = sum(result.*repmat(permute(mac_per_z,[2,3,1]),256,256,1),3);
result_attn_quant = sum(result_quant.*repmat(permute(mac_per_z,[2,3,1]),256,256,1),3);

% a_1_sliced(a_2_sliced > 0.1) = 0;
a_2_sliced(a_2_sliced <= 0.1) = 0;

subplot(2,3,1);
imagesc(a_1_sliced); title('water'); colormap('bone'); colorbar;
subplot(2,3,2);
imagesc(a_2_sliced); title('bone'); colormap('bone'); colorbar;
subplot(2,3,3);
imagesc(result_attn); title('atten seperated'); colormap('bone'); colorbar;
subplot(2,3,4);
imagesc(water_density_sliced+bone_density_sliced); title('density'); colormap('bone'); colorbar;
subplot(2,3,5);
imagesc(mu_mat(1)*a_1_sliced+mu_mat(2)*a_2_sliced); title('attenuation'); colormap('bone'); colorbar;
subplot(2,3,6);
imagesc(result_attn_quant); title('atten seperated quantized'); colormap('bone'); colorbar;
