clear; clc;

% plot_dates = {'220106','220113','220118','220201'};
% photons_label = {'1e8 Photons','5e8 Photons','10e8 Photons','15e8 Photons'};
plot_dates = {'220113','220118','220201','220208'};
photons_label = {'5e8 Photons','10e8 Photons','15e8 Photons','High Res'};
num_loops = length(plot_dates);

plot_images = 1;
plot_scatters = 1;
plot_polar = 1;

images_figure = 1;
norm_images_figure = 2;
scatters_figure = 3;
polar_figure = 4;

for jj=1:num_loops
    
    load([plot_dates{jj},'/reconstruction_atten_anti_grid_Enabled.mat']);
    
    slice = 5;
    figure(images_figure);
    subplot(num_loops,6,1+(jj-1)*6); imagesc(reconstruct_bin_0(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone'); title(photons_label{jj});
    subplot(num_loops,6,2+(jj-1)*6); imagesc(reconstruct_bin_1(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,3+(jj-1)*6); imagesc(reconstruct_bin_2(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,4+(jj-1)*6); imagesc(reconstruct_bin_3(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,5+(jj-1)*6); imagesc(reconstruct_bin_4(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,6+(jj-1)*6); imagesc(reconstruct_bin_5(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone');

    slice = 5;
    figure(norm_images_figure);
    subplot(num_loops,6,1+(jj-1)*6); imagesc(reconstruct_bin_0(:,:,slice)); pbaspect([1,1,1]); colormap('bone'); title(photons_label{jj});
    subplot(num_loops,6,2+(jj-1)*6); imagesc(reconstruct_bin_1(:,:,slice)); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,3+(jj-1)*6); imagesc(reconstruct_bin_2(:,:,slice)); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,4+(jj-1)*6); imagesc(reconstruct_bin_3(:,:,slice)); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,5+(jj-1)*6); imagesc(reconstruct_bin_4(:,:,slice)); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,6+(jj-1)*6); imagesc(reconstruct_bin_5(:,:,slice)); pbaspect([1,1,1]); colormap('bone');

    load([plot_dates{jj},'/GT_I_for_cross_sections.mat']);

    center_idx = 169;
    side_lobe = 60;
    I_GT = permute(I_GT,[2,3,1]);
    I_GT(:,:,3) = I_GT(:,:,3)-I_GT(:,:,1);
    I_GT(1:(center_idx-side_lobe),:,1) = 0;
    I_GT((center_idx+side_lobe):end,:,1) = 0;
    I0_GT = permute(I0_GT,[2,3,1]);
    %%
    detector_types = { ...
    'No Scatter', ...
    'Single Scatter', ...
    'Multi Scatter', ...
    'Single Scatter Compton', ...
    'Single Scatter Rayleigh' ...
    };

    selected_detectors = [2,3,4,5];
    I_GT_selected = I_GT(:,:,selected_detectors);

    sum_energies_selected = squeeze(sum(I_GT_selected,[1,2]));

    num_to_plot = length(selected_detectors);
    if plot_scatters
        figure(scatters_figure);
        for ii=1:num_to_plot
        %     subplot(num_to_plot,2,1+(ii-1)*2);
            subplot(num_to_plot,num_loops,(ii-1)*num_loops+jj);
            very_wide_image = circshift(I_GT_selected(:,:,ii)',size(I_GT,1)/2-center_idx,2);
%             imshow((very_wide_image(:,1:2:end)+very_wide_image(:,2:2:end)),[]);
            imshow((very_wide_image(:,1:4:end)+very_wide_image(:,2:4:end)+very_wide_image(:,3:4:end)+very_wide_image(:,4:4:end)),[]);
            title([detector_types{selected_detectors(ii)},sprintf(' Total Energy: %.2d',sum_energies_selected(ii))])
            axis on;
%             xticks([1 size(very_wide_image,2)/4 size(very_wide_image,2)/2])
            xticks([1 round(size(very_wide_image,2)/8) size(very_wide_image,2)/4])
            xticklabels({'\theta = -\pi','\theta = 0','\theta = \pi'})
            yticks([0])
            yticklabels({''})
        end
    end

    %%
    direct_mask = I_GT(:,:,1) > 0.60*mean(reshape(I0_GT(50:200,10:70,1),1,[]));
    % imshow(direct_mask')

    %%

    % I_GT_selected_minus_mask = I_GT_selected.*(1-repmat(direct_mask,1,1,num_to_plot));
    I_GT_selected_minus_mask = I_GT_selected;

    sums_per_det_per_angle = squeeze(sum(I_GT_selected_minus_mask(:,:,1:num_to_plot),2));

    polar_offset = center_idx*2*pi/size(sums_per_det_per_angle,1);
    polar_vec = mod(linspace(0,2*pi,size(sums_per_det_per_angle,1))-polar_offset,2*pi);
    polar_mat = repmat(polar_vec',1,num_to_plot);

    if plot_polar
        figure(polar_figure)
        subplot(2,num_loops,jj)
        polarplot(polar_mat,(sums_per_det_per_angle));
        if jj==1
            legend(detector_types(selected_detectors),'Location','south')
        end
        title(['Fixed R-Lim - ',photons_label{jj}]);
        rlim([0 12e-5]);
        subplot(2,num_loops,jj+num_loops)
        polarplot(polar_mat,(sums_per_det_per_angle));
%         legend(detector_types(selected_detectors),'Location','south')
        title(['Dynamic R-Lim - ',photons_label{jj}]);
%         semilogy(polar_mat-pi,(sums_per_det_per_angle),'.');
%         legend(detector_types(selected_detectors))
    end

end

%%

selected_detectors = [2,3,4,5];
I_GT_new = I_GT;
I_GT_new(:,:,3) = I_GT(:,:,3)-I_GT(:,:,1);
I_GT_selected = I_GT_new(:,:,selected_detectors);
sums_per_det_per_angle = squeeze(sum(I_GT_selected(:,:,1:num_to_plot),2));

plot_angles = polar_mat(170+(0:753),[3,4,1,2]);
plot_data = sums_per_det_per_angle(170+(0:753),[3,4,1,2]);
sums_per_det_per_angle_no_scat = squeeze(sum(I_GT(:,:,1),2));
plot_data(:,5) = sums_per_det_per_angle_no_scat(170+(0:753));
plot_angles(:,5) = plot_angles(:,4);
clc;
% sums_per_det_per_angle_g = squeeze(sum(I_GT,2));
for angle_to_sum=[30,60,90,120,150,180]
    max_idx = round((angle_to_sum/180)*size(plot_angles,1));
    fprintf('%d & %d\\%% & %d\\%% & %d\\%% & %d\\%% & & %d\\%% & %d\\%% & %d\\%% & %d\\%%\\\\ \\hline \n',angle_to_sum,round(100*sum(plot_data(1:max_idx,1:4))./sum(plot_data(:,1:4))),round(100*sum(plot_data(1:max_idx,1:4))./sum(plot_data(1:max_idx,5))));
end
%%
% Create a polaraxes
ax = polaraxes;
% figure;
polarplot(ax,plot_angles(:,1:4),plot_data(:,1:4));
% polarplot(ax,plot_angles(1:max_idx,:),plot_data(1:max_idx,:));
legend(detector_types([4,5,2,3]),'Location','north')
% set(gca,'rticklabel',[]);
% Update the theta limits of the polaraxes to show just the first quadrant
ax.ThetaLim = [0 180];

%%
energy_csv_d = csvread('220614/all_csv/run_0outputEnergyDep_0.csv',1,0);
count_csv_d = csvread('220614/all_csv/run_0outputCounter_0.csv',1,0);

figure;
plot(sum(energy_csv_d,2)); hold on;
plot(sum(count_csv_d,2));
energy_csv_sum_ns = sum(sum(energy_csv_d,2));
count_csv_sum_ns = sum(sum(count_csv_d,2));

% energy_csv = csvread('220614/all_csv/run_0outputEnergyDep_1.csv',1,0);%-energy_csv;
% count_csv = csvread('220614/all_csv/run_0outputCounter_1.csv',1,0);%-count_csv;
energy_csv = csvread('220614/all_csv/run_0outputEnergyDep_2.csv',1,0)-energy_csv_d;
count_csv = csvread('220614/all_csv/run_0outputCounter_2.csv',1,0)-count_csv_d;

plot(sum(energy_csv,2));
plot(sum(count_csv,2)); hold off;
legend({'E No scat','C No scat','E Single Scat','C Single Scat'})
energy_csv_sum_os = sum(sum(energy_csv,2));
count_csv_sum_os = sum(sum(count_csv,2));

fprintf('Scatterd - E: %d%%. C: %d%%.\n',round(100*energy_csv_sum_os/(energy_csv_sum_os+energy_csv_sum_ns)),round(100*count_csv_sum_os/(count_csv_sum_os+count_csv_sum_ns)));


%%
energy_csv_mat = zeros(1800,80,3);
count_csv_mat = zeros(1800,80,3);
energy_csv_mat(:,:,1) = csvread('220614/all_csv/run_0outputEnergyDep_0.csv',1,0);
count_csv_mat(:,:,1)  = csvread('220614/all_csv/run_0outputCounter_0.csv',1,0);
energy_csv_mat(:,:,2) = csvread('220614/all_csv/run_0outputEnergyDep_1.csv',1,0);
count_csv_mat(:,:,2)  = csvread('220614/all_csv/run_0outputCounter_1.csv',1,0);
energy_csv_mat(:,:,3) = csvread('220614/all_csv/run_0outputEnergyDep_2.csv',1,0);
count_csv_mat(:,:,3)  = csvread('220614/all_csv/run_0outputCounter_2.csv',1,0);

figure;
subplot(6,1,1); imshow(energy_csv_mat(:,:,1)',[]);
subplot(6,1,2); imshow(energy_csv_mat(:,:,2)',[]);
subplot(6,1,3); imshow(energy_csv_mat(:,:,3)',[]);
subplot(6,1,4); imshow(count_csv_mat(:,:,1)',[]);
subplot(6,1,5); imshow(count_csv_mat(:,:,2)',[]);
subplot(6,1,6); imshow(count_csv_mat(:,:,3)',[]);

%%
figure;
subplot(3,1,1); imshow(I_GT(:,:,1)',[]);
subplot(3,1,2); imshow(I_GT(:,:,2)',[]);
subplot(3,1,3); imshow(I_GT(:,:,3)',[]);

%%

% z_vec=[1,6,7,8,15,20];
% e_vec=(0.01:0.00001:0.12)';
% 
% z_mat = repmat(z_vec,length(e_vec),1);
% e_mat = repmat(e_vec,1,length(z_vec));
% 
% res_mat = 2*asin(0.0133*(z_mat.^(1/3))./e_mat)*180/pi;
% res_mat(imag(res_mat)~=0)=180;
% 
% figure;
% plot(e_mat*1000,res_mat)
% title('75%% of all rayleigh scattered rays are under this \theta');
% legend({'H','C','N','O','P','Ca'})
% xlabel('Energy [kEv]');
% ylabel('\theta [^o]');
clear; clc;

I_GT_1 = load('C:\Users\gchem\Geant4\XRayToF_Git\Parsing\Knee_Low_Res\220630\GT_solid_3.mat');
I_GT_2 = load('C:\Users\gchem\Geant4\XRayToF_Git\Parsing\Knee_Low_Res\220630\GT_solid_2.mat');

I_GT_1 = I_GT_1.I_GT;
I_GT_2 = I_GT_2.I_GT;
I_GT_SUM_1 = sum(I_GT_1,3);
I_GT_SUM_2 = sum(I_GT_2,3);

% figure;
% subplot(2,1,1);
% imshow(squeeze(I_GT_1(1,:,:))',[]);
% subplot(2,1,2);
% imshow(squeeze(I_GT_2(1,:,:))',[]);

% to_plot = [4,5];
% figure;
% for ii=1:2
%     ax(ii) = polaraxes;
%     subplot(1,2,ii);
%     polarplot(I_GT_SUM_1(to_plot(ii),:)); hold on;
%     polarplot(I_GT_SUM_2(to_plot(ii),:));
%     ax(ii).ThetaLim = [0 180];
% end

ax_1 = subplot(1,2,1,polaraxes);
polarplot(I_GT_SUM_1(4,:)); hold on;
polarplot(I_GT_SUM_2(4,:));
legend({'Ca','O + Cu'});
% ax_1.ThetaLim = [0 180];
ax_2 = subplot(1,2,2,polaraxes);
polarplot(I_GT_SUM_1(5,:)); hold on;
polarplot(I_GT_SUM_2(5,:));
legend({'Ca','O + Cu'});
% ax_1.ThetaLim = [0 180];
