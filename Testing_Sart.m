clear; clc;
addpath('AIRTools\');

%%
cut_frames = csvread('source_and_det_2\cut_frames_2.csv');

src_id = cut_frames(1:180,1);
det_id_borders = cut_frames(1:180,[2,3])+1;
det_angle_per_source = zeros(length(src_id),401);
for loop_num=1:length(src_id)
    if diff(det_id_borders(loop_num,:))>0
        det_angle_per_source(loop_num,:) = det_id_borders(loop_num,1):det_id_borders(loop_num,2);
    else
        det_angle_per_source(loop_num,:) = [det_id_borders(loop_num,1):1800,1:det_id_borders(loop_num,2)];
    end
end

det_angle_per_source = det_angle_per_source-1;
det_id_per_source = kron(det_angle_per_source,ones(1,80))*80+repmat(0:79',180,401);
%%
row=zeros(968647350,1);
col=zeros(968647350,1);
val=zeros(968647350,1);
cnt = 1;
src_num = 1;
source_reduction = 2;
num_sources_in_vec = 90; 
sources_vec = 1:90; %sort(randperm(180,num_sources_in_vec));

fprintf('Generating A. Working on source %03d/%03d\n',0,num_sources_in_vec);
for source_id = sources_vec
    fprintf('\b\b\b\b\b\b\b\b%03d/%03d\n',src_num,num_sources_in_vec);
    load(sprintf('sparseA_220510/source_id_%03d.mat',source_id));
    [pixel_idx_vec,voxel_idx_vec,ray_length_vec] = find(SparseA);
    ray_length_vec(abs(ray_length_vec)>1)=0;
    len = length(pixel_idx_vec);
    row(cnt:(cnt+len-1)) = pixel_idx_vec+(src_num-1)*401*80;
    col(cnt:(cnt+len-1)) = voxel_idx_vec;
    val(cnt:(cnt+len-1)) = ray_length_vec;
    cnt = cnt+len;
    src_num = src_num + 1;
end

row(cnt:end) = [];
col(cnt:end) = [];
val(cnt:end) = [];

A = sparse(row,col,val,401*80*num_sources_in_vec,256*256*14);
clear row col val SparseA cnt len pixel_idx_vec voxel_idx_vec ray_length_vec;
%%

if true
    load('220510/dect_ground_truth_runs.mat');
    b_det = zeros(80,401,num_sources_in_vec);
    src_num = 1;
    fprintf('Generating b. Working on source %03d/%03d\n',0,num_sources_in_vec);
    for source_id = (sources_vec-1)
        fprintf('\b\b\b\b\b\b\b\b%03d/%03d\n',src_num,num_sources_in_vec);
        source_id_loc = find(src_id==source_id,1);
        div = squeeze(I_GT(source_id_loc,:,:)./I0_GT(source_id_loc,:,:));
        div(div==0) = 1;
        div(div==Inf) = 1;
        div(isnan(div)) = 1;
        log_div = -log(div)';
        b_det(:,:,src_num) = reshape(log_div(det_id_per_source(source_id_loc,:)+1),80,401);
        src_num = src_num + 1;
    end
    b = b_det(:);

    
    K = 200;

    % [X,info,restart] = sart(A,b,K,0.5*ones(256*256*14,1));
    options.lambda = 1.1;
    options.nonneg = true;
    fprintf('Running Reconstruction..\n');
    [X,info,restart] = sart(A/10,b,K,0*ones(256*256*14,1),options);
    fprintf('\b Done.\n');

    X_reshaped = reshape(X,256,256,14);
    SliderImshow(X_reshaped);
%     save('220510/reconstruction_150s.mat','X_reshaped');
end
%%
% load('C:\Users\gchem\Geant4\XRayToF_Git\Parsing\Knee_Low_Res\220405\xcat_reduced.mat')
% 
% b_gt = A*xcat_density(:);
% b_gt_reshaped = reshape(b_gt,80,401,[]);
% SliderImshow(cat(1,b_gt_reshaped/max(b_gt_reshaped,[],'all'),b_det/max(b_det,[],'all')));
% SliderImshow(cat(1,b_det,b_det));

% %%
% X_reshaped = zeros(256,256,14);
% fprintf('Working on source %03d/%d\n',0,180);
% for source_id = 1:180
%     fprintf('\b\b\b\b\b\b\b\b%03d/%d\n',source_id,180);
%     load(sprintf('sparseA_220503/source_id_%03d.mat',source_id));
% 
%     div = squeeze(I_GT(source_id,:,:)./I0_GT(source_id,:,:));
%     div(div==0) = 1;
%     div(div==Inf) = 1;
%     div(isnan(div)) = 1;
%     log_div = -log(div)';
%     log_div = log_div(det_id_per_source(source_id,:)+1)';
%     K = 70;
% 
%     [X,info,restart] = sart(SparseA,log_div,K);
% 
%     X_reshaped = X_reshaped + reshape(X,256,256,14);
% end
%%
if false
    load('220517/ground_truth_runs_with_bins.mat');
    X_reshaped = zeros(256,256,14,6);
    for bin=1:6
        fprintf('Bin %d. \n',bin);
        b_det = zeros(80,401,num_sources_in_vec);
        src_num = 1;
        fprintf('\bGenerating b. Working on source %03d/%03d\n',0,num_sources_in_vec);
        for source_id = (sources_vec-1)
            fprintf('\b\b\b\b\b\b\b\b%03d/%03d\n',src_num,num_sources_in_vec);
            source_id_loc = find(src_id==source_id,1);
            div = squeeze(I_GT_BINS(source_id_loc,bin,:,:)./I0_GT_BINS(source_id_loc,bin,:,:));
            div(div==0) = 1;
            div(div==Inf) = 1;
            div(isnan(div)) = 1;
            log_div = -log(div)';
            b_det(:,:,src_num) = reshape(log_div(det_id_per_source(source_id_loc,:)+1),80,401);
            src_num = src_num + 1;
        end
        b = b_det(:);


        K = 10;

        % [X,info,restart] = sart(A,b,K,0.5*ones(256*256*14,1));
        options.lambda = 1.1;
        options.nonneg = true;
        fprintf('Running Reconstruction..\n');
        [X,info,restart] = sart(A,b,K,0*ones(256*256*14,1),options);
        X_reshaped(:,:,:,bin) = reshape(X,256,256,14);
        fprintf('\b Done.\n');
    %     SliderImshow(X_reshaped);
    %     fprintf('\b Done. Saving bin %d.\n',bin);
    %     save(sprintf('220517/reconstruct_bin_%d.mat',bin),'X_reshaped');
    end

    save('220517/reconstruction_per_bin_150s.mat','X_reshaped');
%     SliderImshow(permute(reshape(permute(X_reshaped,[2,4,1,3]),256*6,256,[]),[2,1,3]))
end