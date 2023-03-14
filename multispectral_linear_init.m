clear; clc;

load('reconstruction_atten_per_bin1.mat')

%%
for slice = 3 %[1,3,6,10]
    per_bin_image(:,:,1) = reconstruct_bin_0(:,:,slice);
    per_bin_image(:,:,2) = reconstruct_bin_1(:,:,slice);
    per_bin_image(:,:,3) = reconstruct_bin_2(:,:,slice);
    per_bin_image(:,:,4) = reconstruct_bin_3(:,:,slice);
    per_bin_image(:,:,5) = reconstruct_bin_4(:,:,slice);
    per_bin_image(:,:,6) = reconstruct_bin_5(:,:,slice);

    figure;
    for ii=1:6
        subplot(2,3,ii);
%         imshow(squeeze(per_bin_image(:,:,ii)),[0,max(per_bin_image,[],'all')]);
        imshow(squeeze(per_bin_image(:,:,ii)));
    end
end
%%
slice = 5;
filter_size = 2;

per_bin_image_filtered = per_bin_image;
figure;
for ii=1:6
    per_bin_image_filtered(:,:,ii) = conv2(per_bin_image(:,:,ii),ones(filter_size)/(filter_size^2),'same');
    subplot(2,3,ii);
    imshow(squeeze(per_bin_image_filtered(:,:,ii)),[0,max(per_bin_image_filtered,[],'all')]);
end

x_pixel = [];
y_pixel = [];


%%
load('ElementXS.mat')
% %%
% [B,I] = sort(ZVec);
% figure;
% plot(Photon(60,I(1:end-2))); hold on;
% plot(Compton(60,I(1:end-2)));
% plot(Rayleigh(60,I(1:end-2))); hold off;
% figure;
% plot(Photon(20:100,24)); hold on;
% plot(Compton(20:100,24));
% plot(Rayleigh(20:100,24)); hold off;
%%

absorbtion_at_60 = Photon(60,:);
scatter_at_60 = Compton(60,:) + Rayleigh(60,:);

Names_vec = {'H','C','N','O','P','K','Ca'};
wanted_z = [1,6,7,8,15,19,20];
indices = zeros(1,length(wanted_z));
for ii=1:length(wanted_z)
    indices(ii) = find(wanted_z(ii)==ZVec,1);
end


figure;
subplot(2,2,1)
scatter(scatter_at_60(indices),absorbtion_at_60(indices),'filled');
grid on;
xlabel('\sigma Scatter [cm^2]');
ylabel('\sigma Absorption [cm^2]');
title('Cross Sections at E=60[KeV]')
dx = 1e-25; dy = 1e-24; % displacement so the text does not overlay the data points
text(scatter_at_60(indices)+dx, absorbtion_at_60(indices)+dy, Names_vec, 'Fontsize', 10);

subplot(2,2,2)
scatter(Rayleigh(60,indices),Compton(60,indices),'filled');
grid on;
xlabel('\sigma Rayleigh [cm^2]');
ylabel('\sigma Compton [cm^2]');
title('Cross Sections at E=60[KeV]')
dx = -1e-25; dy = 5e-25; % displacement so the text does not overlay the data points
text(Rayleigh(60,indices)+dx,Compton(60,indices)+dy, Names_vec, 'Fontsize', 10);

%%
relevant_slice = 5;

data_mat = zeros(size(reconstruct_bin_0,1),size(reconstruct_bin_0,2),6);
for ii=1:6
    data_mat(:,:,ii) = eval(['reconstruct_bin_',num2str(ii-1),'(:,:,',num2str(relevant_slice),')']);
end

bins_edges = [0,20,37,50,62,100,120];
energy_centers = 0.5*(bins_edges(2:end)-bins_edges(1:end-1))+bins_edges(1:end-1);

relevant_bins = 2:5;
cnt=1;
% figure;
for ii=relevant_bins
    subplot(2,4,4+cnt);
    imshow(data_mat(:,:,ii),[0,2.5079])
    rec_images(:,:,cnt) = data_mat(:,:,ii);
    cnt = cnt + 1;
    title(sprintf('Energy bin center: %.1f [KeV]',energy_centers(ii)))
end

%%

% idx     =   1   2   3   4   5   6   7
Names_vec = {'H','C','N','O','P','K','Ca'};
wanted_z  = [ 1 , 6 , 7 , 8 , 15, 19, 20];

energy_centers = [26.5000   41.5000   50.5000   59.0000   70.0000   98.5000];

selected_indices = [1,4,5,7];

selected_elements_names = Names_vec(selected_indices);
selected_elements_z     = wanted_z(selected_indices);
selected_elements_e     = round(energy_centers(relevant_bins));

Na = 6.02214076e23;
xs_mat = zeros(length(selected_elements_e),length(selected_elements_z));
for ii=1:length(selected_elements_z)
    xs_mat(:,ii) = Photon(selected_elements_e,find(selected_elements_z(ii)==ZVec,1))*Na;
end

rec_images_vecs = reshape(rec_images,[],length(selected_elements_e));


% rec_Ck_vecs = (((xs_mat)^-1)*(rec_images_vecs.')).';
rec_Ck_vecs = (((xs_mat+eye(4)*1e-2)^-1)*(rec_images_vecs.')).';
rec_Ck = reshape(rec_Ck_vecs,size(rec_images,1),size(rec_images,2),[]);

figure;
for ii=1:size(rec_Ck,3)
    subplot(2,2,ii);
    imshow(rec_Ck(:,:,ii),[])
    title(sprintf('Material: %s',selected_elements_names{ii}))
end
