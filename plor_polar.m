clear;clc;

load('220201\reconstruction_atten_anti_grid_Enabled.mat')

slice_to_plot = 5;
per_bin_image(:,:,1) = reconstruct_bin_0(:,:,slice_to_plot);
per_bin_image(:,:,2) = reconstruct_bin_1(:,:,slice_to_plot);
per_bin_image(:,:,3) = reconstruct_bin_2(:,:,slice_to_plot);
per_bin_image(:,:,4) = reconstruct_bin_3(:,:,slice_to_plot);
per_bin_image(:,:,5) = reconstruct_bin_4(:,:,slice_to_plot);
per_bin_image(:,:,6) = reconstruct_bin_5(:,:,slice_to_plot);
%%
x_mat = repmat((1:60)-30.5,60,1);
y_mat = repmat((1:60)'-30.5,1,60);
[theta,rho] = cart2pol(x_mat,y_mat);
rho_vec = unique(rho)';
theta_vec = unique(theta)';

num_rs = 83;
num_thetas = 360;
old_image = reconstruct_with_asg(:,:,5);
new_image = ImToPolar (old_image, 0, 1, num_rs, num_thetas)';

figure;
subplot(1,2,1); imagesc(old_image); pbaspect([1,1,1]); colormap('bone')
subplot(2,2,2); imagesc(new_image);
subplot(2,2,4); plot(sum(new_image)); hold on;
ylim(lines_y(:,1));
xlim([1,num_rs]);
lines_x = repmat(find(islocalmin(sum(new_image))),2,1);
lines_y = repmat([min(sum(new_image));max(sum(new_image))],1,size(lines_x,2));
plot(lines_x,lines_y,'k');
lines_x = repmat(find(islocalmax(sum(new_image))),2,1);
lines_y = repmat([min(sum(new_image));max(sum(new_image))],1,size(lines_x,2));
plot(lines_x,lines_y,'r');
legend({'Sum along R','Min','Max'});
%%
figure(5)
for ii=1:6
    subplot(2,3,ii);
%     imagesc(per_bin_image(:,:,ii),[0,1.3]); colormap('bone');
    imagesc(per_bin_image(:,:,ii)); colormap('bone');
    axis off;
    title(sprintf('Bin %d',ii));
end