clear; clc;
load('I0_for_gt_and_for_forward.mat');


%%
for ii=1:15
    subplot(15,2,ii*2-1);
    imagesc(squeeze(I0_GT(ii,:,:))')
    subplot(15,2,ii*2);
    imagesc(squeeze(I0_forward(ii,:,:))')
end
%%
I0_GT_relevant = I0_GT(1:15,:,:);
I0_GT_relevant(I0_GT_relevant == 0) = inf;

size_short = size(I0_forward,3);
size_long = size(I0_forward,2);

I0_forward_piled = zeros(size_short*15,size_long);
I0_GT_piled = zeros(size_short*15,size_long);
I0_GT_relevant_piled = zeros(size_short*15,size_long);

for ii=1:15
    I0_forward_piled((1:size_short)+size_short*(ii-1),:) = squeeze(I0_forward(ii,:,:))';
    I0_GT_piled((1:size_short)+size_short*(ii-1),:) = squeeze(I0_GT(ii,:,:))';
    I0_GT_relevant_piled((1:size_short)+size_short*(ii-1),:) = squeeze(I0_GT_relevant(ii,:,:))';
end

j_ratio = I0_forward_piled./I0_GT_relevant_piled;

subplot(1,3,1);
imagesc(I0_GT_piled)
subplot(1,3,2);
imagesc(I0_forward_piled)
subplot(1,3,3);
imagesc(j_ratio)

%%
gt_avg_vec = mean(I0_GT(1:15,:,:),[2,3])';
fw_avg_vec = mean(I0_forward(1:15,:,:),[2,3])';

%%

small_j_ratio = medfilt2(j_ratio(1:80,15:300));
imagesc(small_j_ratio)
