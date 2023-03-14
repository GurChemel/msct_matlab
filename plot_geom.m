clear;clc;

src_loc = csvread('source_and_det_2\src_loc.csv');
src_orient = csvread('source_and_det_2\src_orient.csv');
det_loc = csvread('source_and_det_2\det_loc.csv');
cut_frames = csvread('source_and_det_2\cut_frames_2.csv');
%%
det_to_plot = 40;
scatter3(src_loc(:,1),src_loc(:,2),src_loc(:,3),'.'); hold on;
scatter3(det_loc(1:det_to_plot,1),det_loc(1:det_to_plot,2),det_loc(1:det_to_plot,3),'.');


%%
figure;
src_id = [1,20, 90];
det_id = (src_id-1)*10+1;
mag = 265+267;
det_loc_flat = det_loc(1:80:end,:);
scatter(src_loc(src_id,1),src_loc(src_id,2),'.'); hold on;
scatter([0;det_loc_flat(det_id,1)],[0;det_loc_flat(det_id,2)],'x'); hold on;
% quiver(src_loc(src_id,1),src_loc(src_id,2),mag*src_orient(src_id,1),mag*src_orient(src_id,3),0); hold on;
quiver(src_loc(src_id,1),src_loc(src_id,2),mag*src_orient(src_id,4),mag*src_orient(src_id,6),0); hold on;
% quiver(src_loc(src_id,1),src_loc(src_id,2),mag*src_orient(src_id,1)+mag*src_orient(src_id,4),mag*src_orient(src_id,3)+mag*src_orient(src_id,6),0); hold on;

%%
figure;
num_plots = 6;
rand_indices = randi(180,1,num_plots);
src_id = cut_frames(rand_indices,1);
det_id = cut_frames(rand_indices,[2,3]);
mag = 265+267;
det_loc_flat = det_loc(1:80:end,:);
for loop_num=1:num_plots
    subplot(2,3,loop_num);
    % src
    scatter(src_loc(src_id(loop_num),1),src_loc(src_id(loop_num),2),'.'); hold on;
    quiver(src_loc(src_id(loop_num),1),src_loc(src_id(loop_num),2),mag*src_orient(src_id(loop_num),4),mag*src_orient(src_id(loop_num),6),0); hold on;
    % det
    if diff(det_id(loop_num,:))>0
        plot_det_id = det_id(loop_num,1):det_id(loop_num,2);
    else
        plot_det_id = [det_id(loop_num,1):1800,1:det_id(loop_num,2)];
    end
    scatter([0;det_loc_flat(:,1)],[0;det_loc_flat(:,2)],'.'); pbaspect([1,1,1]); hold on;
    scatter([0;det_loc_flat(plot_det_id,1)],[0;det_loc_flat(plot_det_id,2)],'.'); pbaspect([1,1,1]); hold off;
end
%%

plot_names = {'cone','cone_vec','cone_vec_fixed'};
num_loops = length(plot_names);

for jj=1:num_loops
    
    load(['220208_3\reconstruction_atten_anti_grid_Enabled_',plot_names{jj},'.mat']);
    title_vec = plot_names{jj}; title_vec(title_vec=='_') = ' ';

    slice = 25;
    figure(2);
    subplot(1,num_loops,jj); imagesc(reconstruct(:,:,slice)); pbaspect([1,1,1]); colormap('bone'); title(title_vec);
    
    figure(3);
    subplot(num_loops,6,1+(jj-1)*6); imagesc(reconstruct_bin_0(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone'); title(title_vec);
    subplot(num_loops,6,2+(jj-1)*6); imagesc(reconstruct_bin_1(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,3+(jj-1)*6); imagesc(reconstruct_bin_2(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,4+(jj-1)*6); imagesc(reconstruct_bin_3(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,5+(jj-1)*6); imagesc(reconstruct_bin_4(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,6+(jj-1)*6); imagesc(reconstruct_bin_5(:,:,slice),[0,1.3]); pbaspect([1,1,1]); colormap('bone');

    figure(4);
    subplot(num_loops,6,1+(jj-1)*6); imagesc(reconstruct_bin_0(:,:,slice)); pbaspect([1,1,1]); colormap('bone'); title(title_vec);
    subplot(num_loops,6,2+(jj-1)*6); imagesc(reconstruct_bin_1(:,:,slice)); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,3+(jj-1)*6); imagesc(reconstruct_bin_2(:,:,slice)); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,4+(jj-1)*6); imagesc(reconstruct_bin_3(:,:,slice)); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,5+(jj-1)*6); imagesc(reconstruct_bin_4(:,:,slice)); pbaspect([1,1,1]); colormap('bone');
    subplot(num_loops,6,6+(jj-1)*6); imagesc(reconstruct_bin_5(:,:,slice)); pbaspect([1,1,1]); colormap('bone');
end