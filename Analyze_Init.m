clear; clc;
% load('220918/reconstruction_atten_per_bin0.mat');
% load('220918/reconstruction_full_spectrum.mat');
% load('220918/ground_truth_runs_with_bins.mat');
% slice = 5;

% load('221025/reconstruction_atten_per_bin0.mat');
% load('221025/reconstruction_full_spectrum.mat');
% load('221025/ground_truth_run_with_bins.mat');
% slice = 5;

% load('221110/reconstruction_atten_per_bin0.mat');
% load('221110/reconstruction_full_spectrum.mat');
% load('221110/ground_truth_run_with_bins.mat');
% slice = 1;

% load('221117/reconstruction_atten_per_bin0.mat');
% load('221117/reconstruction_full_spectrum.mat');
% load('221117/ground_truth_run_with_bins.mat');
% slice = 5;

load('230202/reconstruction_full_spectrum.mat');
slice = 1;

%%
figure;
subplot(3,1,1);
imagesc(I_GT_single');
colormap('bone');
subplot(3,1,2);
imagesc(I_GT_scatter');
colormap('bone');
subplot(3,1,3);
imagesc(I_GT_asg');
colormap('bone');
%%
I_GT_direct_norm = (I_GT_single/max(I_GT_single,[],'all')).^0.5;
I_GT_scatter_norm = (I_GT_scatter/max(I_GT_scatter,[],'all')).^0.5;
I_GT_asg_norm = (I_GT_asg/max(I_GT_asg,[],'all')).^0.5;

figure;
subplot(3,1,1);
imshow(I_GT_direct_norm');
% colormap('bone');
subplot(3,1,2);
imshow(I_GT_scatter_norm');
% colormap('bone');
subplot(3,1,3);
imshow(I_GT_asg_norm');
% colormap('bone');

imshow(I_GT_asg_norm' - I_GT_direct_norm',[]);

%%
filename = "I_GT_Animated.gif"; % Specify the output file name
I_GT_norm = (I_GT/max(I_GT,[],'all')).^0.5; % Gamma correction
I_GT_norm = uint8(255*I_GT_norm);
for idx = 1:6:180
    if idx == 1
        imwrite(squeeze(I_GT_norm(idx,:,:)).',filename,"gif","LoopCount",Inf,"DelayTime",0.1);
    else
        imwrite(squeeze(I_GT_norm(idx,:,:)).',filename,"gif","WriteMode","append","DelayTime",0.1);
    end
end
imwrite(squeeze(I_GT_norm(1,:,:)).',"I_GT_First.png");
%%
imwrite(squeeze(I_GT_norm(1,:,:)).',"ImagesForThesis\I_GT_First.png");
imwrite(squeeze(I_GT_norm(46,:,:)).',"ImagesForThesis\I_GT_45th.png");

%%


figure;
imshow(cat(2,reconstruct_no_asg(:,:,slice),reconstruct_with_asg(:,:,slice)),[0,0.51]); title('Without ASG | With ASG')
%%
reconstruct_no_asg_norm = (reconstruct_no_asg/max(reconstruct_no_asg,[],'all')).^0.5;
reconstruct_with_asg_norm = (reconstruct_with_asg/max(reconstruct_with_asg,[],'all')).^0.5;

figure;
subplot(1,2,1);
imshow(reconstruct_no_asg_norm(:,:,slice)); title('without ASG')
% colormap('bone');
subplot(1,2,2);
imshow(reconstruct_with_asg_norm(:,:,slice)); title('with ASG')
% colormap('bone');

% % xcat_density_norm = (xcat_density/max(xcat_density,[],'all')).^0.5;
% % subplot(1,2,2);
% % imshow(xcat_density_norm(:,:,slice)); title('with ASG')
% 
% gt_linear_atten_full_norm = (gt_linear_atten_full/max(gt_linear_atten_full,[],'all')).^0.5;
% subplot(1,2,2);
% imshow(gt_linear_atten_full_norm(:,:,slice)); title('with ASG')
%%
imwrite(reconstruct_no_asg_norm(:,:,slice),"ImagesForPresentations/recon_noised.png");
imwrite(gt_linear_atten_full_norm(:,:,slice),"ImagesForPresentations/recon_clean.png");

%%

subplot(1,2,2);
% figure;
imshow(reconstruct_with_asg_norm(:,:,slice)); title('with ASG')
% colormap('bone');
%%
figure;
for bin=1:6
    image_values(:,:,bin) = eval(sprintf('reconstruct_bin_%d(:,:,%d)',bin-1,slice));
end
image_values = image_values/max(image_values,[],'all');
image_values = (image_values.^0.5);
for bin=1:6
    subplot(2,3,bin);
    imshow(image_values(:,:,bin));
%     imwrite(image_values(:,:,bin),sprintf('ImagesForPresentations/gamma_corrected_bin_%d.jpg',bin));
%     imagesc(image_values(:,:,bin));
%     colormap('bone');
end
%%
filename = "gamma_corrected_bins_Animated.gif"; % Specify the output file name
image_values_uint8 = uint8(255*image_values);
for bin = 1:6
    if bin == 1
        imwrite(squeeze(image_values_uint8(:,:,bin)),filename,"gif","LoopCount",Inf,"DelayTime",0.5);
    else
        imwrite(squeeze(image_values_uint8(:,:,bin)),filename,"gif","WriteMode","append","DelayTime",0.5);
    end
end

%%
figure;
for bin=1:6
    subplot(2,3,bin);
    imagesc(squeeze(I_GT_BINS(bin,:,:))');
    colormap('bone');
end
%%
labels = zeros(size(reconstruct_with_asg));
for slice_num=1:size(reconstruct_with_asg,3)
    input_image = zeros(150,150,6);
    for bin=1:6
        input_image(:,:,bin) = eval(sprintf('reconstruct_bin_%d(:,:,%d)',bin-1,slice_num));
    end
    smothness_factor = 0.02;
    labels(:,:,slice_num) = GraphCutPerBinWrapper(input_image, smothness_factor, 1);
end
figure;
subplot(1,2,1); imshow(reconstruct_with_asg(:,:,slice),[]);
subplot(1,2,2); imagesc(labels(:,:,slice)); pbaspect([1,1,1]);

bone_recon_with_asg = reconstruct_with_asg.*(labels==2);
water_recon_with_asg = reconstruct_with_asg.*(labels==3);

SliderImshow(cat(2,bone_recon_with_asg,water_recon_with_asg))

%%

h = volshow(bone_recon_with_asg,config);
viewer = h.Parent;
hFig = viewer.Parent;
drawnow

I = getframe(hFig);
[indI,cm] = rgb2ind(I.cdata,256);
imwrite(indI,cm,'tansparent_test.png', 'png', 'transparency', [0.3020    0.7490    0.9294]) %,"gif",Loopcount=inf,DelayTime=0)

% imwrite(bitmapData, 'a.png', 'png', 'transparency', backgroundColor)

%     % Write the frame to the GIF file
%     if idx==1
%         % Do nothing. The first frame displays only the viewer, not the
%         % volume.
%     elseif idx == 2
%         imwrite(indI,cm,filename,"gif",Loopcount=inf,DelayTime=0)
%     else
%         imwrite(indI,cm,filename,"gif",WriteMode="append",DelayTime=0)
%     end


%%
% figure;
% imagesc(labels(:,:,slice))
figure;
imagesc(bone_recon_with_asg(:,:,slice)); colormap('bone');
figure;
imagesc(water_recon_with_asg(:,:,slice)); colormap('bone');


%%
clear; clc;
load('220920\xcat_reduced.mat')
load('220920\Single_0init_vol.mat')

slice = 35;
rec_density = init_vol_density(:,:,slice);
xct_density = xcat_density(:,:,slice);

max_val = max(cat(3,rec_density,xct_density),[],'all');
min_val = min(cat(3,rec_density,xct_density),[],'all');

figure;
subplot(2,2,1);
imshow(rec_density,[min_val max_val]);
subplot(2,2,2);
imshow(xct_density,[min_val max_val]);
subplot(2,2,3);
imshow((xct_density-rec_density),[]);
subplot(2,2,4);
imshow(abs(xct_density-rec_density),[]);

%%
clear;clc;

for phantom_id = 0:1
    if phantom_id==0
        load('221208\xcat_reduced.mat')
        load('221208\ground_truth_run_with_bins.mat')
        num_photons = 10000000; % *10
        run_factor = 180;
    else
        load('230202\xcat_reduced.mat')
        load('230202\ground_truth_runs_no_bins.mat')
        I_GT_scatter = I_GT;
        I0_GT_scatter = I0_GT;
        num_photons = 10000000; %*100
        run_factor = 1;
    end
    % xcat_density [g/cm^3].
    % 0.001*xcat_density [kg/cm^3].
    % Voxel size: 1x1x1 mm = 0.1*0.1*0.1 cm
    % 0.001*xcat_density*0.001 [kg].
    sum_kg = sum(0.001*xcat_density*0.001,'all');
    disp(['Weight [Kg] is: ',num2str(sum_kg)])

    I_GT_REAL=I_GT_scatter*run_factor*num_photons;
    I0_GT_REAL=I0_GT_scatter*run_factor*num_photons;

    % 1 kEv = 1.60218e-16 Joule
    absorbed_energy_kev = (sum(I0_GT_REAL,'all')-sum(I_GT_REAL,'all'));
    disp(['Absorbed Energky [KeV] is: ',num2str(absorbed_energy_kev)])
    absorbed_energy_joule = absorbed_energy_kev*(1.60218e-16);

    % Absorbed dose: 1 Gy = 1000 mSv = 1 J/Kg
    % 1 mSv = 0.001 J/Kg
    % 1 J/Kg = 1000 mSv
    average_mSv = 1000*absorbed_energy_joule/sum_kg;
    disp(['Effective dose is: ',num2str(average_mSv)])
end


%%
clear;clc;
load('params_for_error_calc.mat')
Materials = load('Materials.mat');

load('221208\init_vol.mat')
load('221208\xcat_reduced.mat')
load('221208\reconstruction_full_spectrum.mat')
load('221208\reconstruction_atten_per_bin0.mat')
slice = 13;

% load('230204\init_vol.mat')
% load('230202\xcat_reduced.mat')
% load('230204\reconstruction_full_spectrum.mat')
% load('230204\reconstruction_atten_per_bin0.mat')
% slice = 1;

%%
figure;
for bin=1:6
    subplot(2,3,bin);
    cur_slice = eval(sprintf('reconstruct_bin_%d(:,:,%d)',bin-1,slice));
    cur_slice = cur_slice./max(cur_slice,[],'all');
%     cur_slice = cur_slice.^0.8;
    imwrite(cur_slice,sprintf('ImagesForThesis/opt_gct_lin_recon_bin_%d.jpg',bin));
    imshow(cur_slice,[]);
end
%%
init_slice = no_asg_init_vol_density(:,:,slice);
xcat_slice = xcat_density(:,:,slice);
init_ids = xcat_id(:,:,slice);

% combined_image = cat(2,,xcat_density(:,:,slice));
% imshow(combined_image,[])
bone_vals = reshape(init_slice(init_ids==38),[],1);
marrow_vals = reshape(init_slice(init_ids==31),[],1);
tissue_vals = reshape(init_slice(init_ids==6),[],1);
blood_vals = reshape(init_slice(init_ids==2),[],1);

figure;
subplot(2,2,1)
plot(bone_vals); hold on;
plot(reshape(xcat_slice(init_ids==38),[],1)); hold on;
subplot(2,2,2)
plot(marrow_vals); hold on;
plot(reshape(xcat_slice(init_ids==31),[],1)); hold on;
subplot(2,2,3)
plot(tissue_vals); hold on;
plot(reshape(xcat_slice(init_ids==6),[],1)); hold on;
subplot(2,2,4)
plot(blood_vals); hold on;
plot(reshape(xcat_slice(init_ids==2),[],1)); hold on;

%%
addpath('GCMex/')

reconstruct_per_bin(:,:,:,1) = reconstruct_bin_0;
reconstruct_per_bin(:,:,:,2) = reconstruct_bin_1;
reconstruct_per_bin(:,:,:,3) = reconstruct_bin_2;
reconstruct_per_bin(:,:,:,4) = reconstruct_bin_3;
reconstruct_per_bin(:,:,:,5) = reconstruct_bin_4;
reconstruct_per_bin(:,:,:,6) = reconstruct_bin_5;

% load('220510\reconstruction_150s.mat')
% reconstruct = X_reshaped*2.8;
% load('220517\reconstruction_per_bin_150s.mat')
% reconstruct_per_bin = X_reshaped*2.8;

LoadGParams;
smothness_factor = 0.02;

materials_recovered_gct = GraphCut3dPerBinWrapper(squeeze(reconstruct_per_bin), smothness_factor, 3);
% materials_recovered_gct = GraphCut3dPerBinWrapper(squeeze(reconstruct_per_bin), smothness_factor, 1);
SliderImagesc(materials_recovered_gct)
%%
Materials_names = {'Air','Water','Muscle','dryribwater','Redmarrow','Iodineblood','Adiposetissue'};
active_materials = [1, 2, 3, 4, 5];
% active_materials = [1, 4, 7];
Materials_names = Materials_names(active_materials);
graphcuts_fractional_mass_full = zeros(size(xcat_fractional_mass));
for ii=1:size(materials_recovered_gct,1)
    for jj=1:size(materials_recovered_gct,2)
        for kk=1:size(materials_recovered_gct,3)
            vec = Materials.(Materials_names{materials_recovered_gct(ii,jj,kk)});
            graphcuts_fractional_mass_full(:,ii,jj,kk) = vec(1:27);
        end
    end
end
graphcuts_density = no_asg_init_vol_density;
graphcuts_fractional_mass = graphcuts_fractional_mass_full(1+[0,7,13,18],:,:,:);
graphcuts_fractional_mass = graphcuts_fractional_mass./(repmat(sum(graphcuts_fractional_mass,1),4,1,1,1));
%%
zzz=1+[0,5,6,7,13,18];
fprintf('               &   H    &   C    &   N    &   O    &   P    &   Ca   \\\\ \\hline\n');
for jj=1:length(Materials_names)
    vec = Materials.(Materials_names{jj});
    vec = vec(zzz)./sum(vec(zzz));
    fprintf('%-14s',Materials_names{jj});
    fprintf(' & %6.4f',vec);
    fprintf(' \\\\\n');
end
% graphcuts_fractional_mass_for_thesis = graphcuts_fractional_mass_full(1+[0,5,6,7,13,18],:,:,:);
% graphcuts_fractional_mass_for_thesis = graphcuts_fractional_mass_for_thesis./(repmat(sum(graphcuts_fractional_mass_for_thesis,1),6,1,1,1));
% fractional_mass_mat = unique(reshape(graphcuts_fractional_mass_for_thesis,6,[])','rows');
%%
element_location_dict_in_xcat   = struct('H', 0, 'C', 5, 'N', 6, 'O', 7, 'P', 13, 'Ca', 18);
element_location_for_error_dict = struct('H', 0, 'C', 1, 'N', 2, 'O', 3, 'P',  4, 'Ca',  5);

asg_linear_att_per_element = zeros(size(asg_init_vol_fractional_mass));
no_asg_linear_att_per_element = zeros(size(asg_init_vol_fractional_mass));
gt_linear_att_per_element = zeros(size(asg_init_vol_fractional_mass));
graphcuts_linear_att_per_element = zeros(size(asg_init_vol_fractional_mass));

active_elements = {'H','O','P','Ca'};
for ind=1:length(active_elements)
    el = active_elements{ind};
    location_in_xcat = element_location_dict_in_xcat.(el)+1; % +1 for matlab
    location_in_mu_mat = element_location_for_error_dict.(el)+1; % +1 for matlab
    
    asg_dens_per_element = squeeze(asg_init_vol_fractional_mass(ind, :, :, :)) .* asg_init_vol_density;
    no_asg_dens_per_element = squeeze(no_asg_init_vol_fractional_mass(ind, :, :, :)) .* no_asg_init_vol_density;
    gt_dens_per_element = squeeze(xcat_fractional_mass(location_in_xcat, :, :, :)) .* xcat_density;
    graphcuts_dens_per_element = squeeze(graphcuts_fractional_mass_full(location_in_xcat, :, :, :)) .* graphcuts_density;
    for e_idx=1:length(source_energy_vec)
        energy_p = source_energy_vec(e_idx);
        atten_cm2g = atten_per_element_cm2g_kev(e_idx, location_in_mu_mat);
        asg_linear_att_per_element(ind,:,:,:) = asg_linear_att_per_element(ind,:,:,:) + permute(energy_p * asg_dens_per_element * atten_cm2g,[4,1,2,3]);
        no_asg_linear_att_per_element(ind,:,:,:) = no_asg_linear_att_per_element(ind,:,:,:) + permute(energy_p * no_asg_dens_per_element * atten_cm2g,[4,1,2,3]);
        gt_linear_att_per_element(ind,:,:,:) = gt_linear_att_per_element(ind,:,:,:) + permute(energy_p * gt_dens_per_element * atten_cm2g,[4,1,2,3]);
        graphcuts_linear_att_per_element(ind,:,:,:) = graphcuts_linear_att_per_element(ind,:,:,:) + permute(energy_p * graphcuts_dens_per_element * atten_cm2g,[4,1,2,3]);
    end
    
end

asg_linear_atten_full = squeeze(sum(asg_linear_att_per_element));
no_asg_linear_atten_full = squeeze(sum(no_asg_linear_att_per_element));
gt_linear_atten_full = squeeze(sum(gt_linear_att_per_element));
graphcuts_linear_atten_full = squeeze(sum(graphcuts_linear_att_per_element));
% %%
% asg_linear_atten_full = asg_linear_atten_full(:,:,slice);
% no_asg_linear_atten_full = no_asg_linear_atten_full(:,:,slice);
% gt_linear_atten_full = gt_linear_atten_full(:,:,slice);
% graphcuts_linear_atten_full = graphcuts_linear_atten_full(:,:,slice);

asg_epsilon = sum(abs(gt_linear_atten_full-asg_linear_atten_full),'all')/sum(abs(gt_linear_atten_full),'all');
asg_delta   = (sum(abs(asg_linear_atten_full),'all')-sum(abs(gt_linear_atten_full),'all'))/sum(abs(gt_linear_atten_full),'all');
no_asg_epsilon = sum(abs(gt_linear_atten_full-no_asg_linear_atten_full),'all')/sum(abs(gt_linear_atten_full),'all');
no_asg_delta   = (sum(abs(no_asg_linear_atten_full),'all')-sum(abs(gt_linear_atten_full),'all'))/sum(abs(gt_linear_atten_full),'all');
graphcuts_epsilon = sum(abs(gt_linear_atten_full-graphcuts_linear_atten_full),'all')/sum(abs(gt_linear_atten_full),'all');
graphcuts_delta   = (sum(abs(graphcuts_linear_atten_full),'all')-sum(abs(gt_linear_atten_full),'all'))/sum(abs(gt_linear_atten_full),'all');

fprintf('With    ASG. Epsilon: %.3f, Delta: %.3f\n',asg_epsilon,asg_delta);
fprintf('Without ASG. Epsilon: %.3f, Delta: %.3f\n',no_asg_epsilon,no_asg_delta);
fprintf('Graph  Cuts. Epsilon: %.3f, Delta: %.3f\n',graphcuts_epsilon,graphcuts_delta);

SliderImshow(cat(2,gt_linear_atten_full,graphcuts_linear_atten_full,no_asg_linear_atten_full))

%%
for sl=[2,13]
    imwrite(kron(gt_linear_atten_full(:,:,sl),ones(2)),sprintf('ImagesForThesis/gct_full_atten_gt_slice_%d.jpg',sl));
    imwrite(kron(graphcuts_linear_atten_full(:,:,sl),ones(2)),sprintf('ImagesForThesis/gct_full_atten_gct_slice_%d.jpg',sl));
    imwrite(kron(no_asg_linear_atten_full(:,:,sl),ones(2)),sprintf('ImagesForThesis/gct_full_atten_naive_slice_%d.jpg',sl));
end
