clear; clc;
load('reconstruction_full_spectrum.mat')

slice = 5;
figure;
subplot(2,3,1);
imagesc(squeeze(reconstruct_no_asg(:,:,slice)));
colormap('bone'); pbaspect([1,1,1]);
title('No ASG')
subplot(2,3,2);
imagesc(squeeze(reconstruct_with_asg(:,:,slice)));
colormap('bone'); pbaspect([1,1,1]);
title('With ASG')
subplot(2,3,3);
imagesc(squeeze(reconstruct_with_asg_2(:,:,slice)));
colormap('bone'); pbaspect([1,1,1]);
title('With ASG 2')
    
subplot(2,1,2);
imagesc(squeeze([reconstruct_no_asg(:,:,slice), reconstruct_with_asg(:,:,slice), reconstruct_with_asg_2(:,:,slice)]));
colormap('bone'); pbaspect([2,1,1]);
title('No | With | With 2 - ASG')
