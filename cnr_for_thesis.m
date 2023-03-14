clear;clc;

load('230204/reconstruction_full_spectrum.mat');
slice = 1;

asg_slice = reconstruct_with_asg(:,:,slice);
no_asg_slice = reconstruct_no_asg(:,:,slice);

x_start = 150;
x_end = x_start+150;
y_start = 170;
y_end = y_start+110;
window_1 = zeros(size(asg_slice));
window_1(x_start:x_end,y_start:y_end)=1;
window_1((x_start+1):(x_end-1),(y_start+1):(y_end-1))=0;

asg_win_1 = asg_slice((x_start+1):(x_end-1),(y_start+1):(y_end-1));
no_asg_win_1 = no_asg_slice((x_start+1):(x_end-1),(y_start+1):(y_end-1));

x_start = 360;
x_end = x_start+140;
y_start = 240;
y_end = y_start+170;
window_2 = zeros(size(asg_slice));
window_2(x_start:x_end,y_start:y_end)=1;
window_2((x_start+1):(x_end-1),(y_start+1):(y_end-1))=0;

asg_win_2 = asg_slice((x_start+1):(x_end-1),(y_start+1):(y_end-1));
no_asg_win_2 = no_asg_slice((x_start+1):(x_end-1),(y_start+1):(y_end-1));

asg_slice=asg_slice + window_1 + window_2;
no_asg_slice=no_asg_slice + window_1 + window_2;


figure;
subplot(1,2,1);
imshow(asg_slice,[0,0.51]); title('With ASG')
subplot(1,2,2);
imshow(no_asg_slice,[0,0.51]); title('Without ASG')

imwrite(asg_slice,"ImagesForThesis\Lin_With_ASG.png");
imwrite(no_asg_slice,"ImagesForThesis\Lin_Without_ASG.png");

window_1_regions = asg_win_1>=0.12;
window_2_regions = asg_win_2>=0.3;

% zzz=asg_win_2;
% zzz(zzz<0.3)=0;
% imshow(zzz,[0,0.51]); title('ASG Window 1');

cnr_asg_win_1 = (mean(asg_win_1(window_1_regions))-mean(asg_win_1(window_1_regions==0)))/(sqrt(((std(asg_win_1(window_1_regions))^2)+(std(asg_win_1(window_1_regions==0))^2))*0.5));
cnr_asg_win_2 = (mean(asg_win_2(window_2_regions))-mean(asg_win_2(window_2_regions==0)))/(sqrt(((std(asg_win_2(window_2_regions))^2)+(std(asg_win_2(window_2_regions==0))^2))*0.5));
cnr_no_asg_win_1 = (mean(no_asg_win_1(window_1_regions))-mean(no_asg_win_1(window_1_regions==0)))/(sqrt(((std(no_asg_win_1(window_1_regions))^2)+(std(no_asg_win_1(window_1_regions==0))^2))*0.5));
cnr_no_asg_win_2 = (mean(no_asg_win_2(window_2_regions))-mean(no_asg_win_2(window_2_regions==0)))/(sqrt(((std(no_asg_win_2(window_2_regions))^2)+(std(no_asg_win_2(window_2_regions==0))^2))*0.5));

figure;
subplot(2,2,1);
imshow(asg_win_1,[0,0.51]); title(['ASG Window 1. CNR: ',num2str(cnr_asg_win_1)]);
subplot(2,2,2);
imshow(no_asg_win_1,[0,0.51]); title(['NO ASG Window 1. CNR: ',num2str(cnr_no_asg_win_1)]);
subplot(2,2,3);
imshow(asg_win_2,[0,0.51]); title(['ASG Window 2. CNR: ',num2str(cnr_asg_win_2)]);
subplot(2,2,4);
imshow(no_asg_win_2,[0,0.51]); title(['NO ASG Window 2. CNR: ',num2str(cnr_no_asg_win_2)]);


