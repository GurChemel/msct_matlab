clear; clc;
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 80);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Built Database:
num_detectors = 21;
Output_counters = zeros(80,1508,20);

% Import the data and convert:
% run0outputCounter0 = readtable("210923\run_0outputCounter_0.csv", opts);
% run0outputCounter0 = table2array(run0outputCounter0)';
% run0outputCounter2 = readtable("210923\run_0outputCounter_2.csv", opts);
% run0outputCounter2 = table2array(run0outputCounter2)';
for ii=1:num_detectors
    Output_counters(:,:,ii) = table2array(readtable(sprintf('211001\\run_0outputCounter_%d.csv',ii-1), opts))';
end

% Clear temporary variables
clear opts

%% Start calculating

subplot(3,1,1);
imagesc(Output_counters(:,:,1)); colormap('bone');
subplot(3,1,2);
imagesc(Output_counters(:,:,3)); colormap('bone');
subplot(3,1,3);
imagesc(Output_counters(:,:,3)-Output_counters(:,:,1)); colormap('bone');


x_vec = 30:280;

direct_image = Output_counters(:,x_vec,1)/max(Output_counters(:,x_vec,3),[],'all');
direct_and_global_image = Output_counters(:,x_vec,3)/max(Output_counters(:,x_vec,3),[],'all');
global_image = direct_and_global_image-direct_image;

min_val = min([direct_image,direct_and_global_image,global_image],[],'all');
max_val = max([direct_image,direct_and_global_image,global_image],[],'all');

direct_image_norm = imadjust(direct_image,[min_val, max_val],[1,0]);
direct_and_global_image_norm = imadjust(direct_and_global_image,[min_val, max_val],[1,0]);
global_image_norm = direct_image_norm-direct_and_global_image_norm;

write_im_for_pres(direct_image_norm,'Counting Only Direct Photons',2)
write_im_for_pres(direct_and_global_image_norm,'Counting Direct and Global Photons',2)
write_im_for_pres(global_image_norm,'Counting Only Global Photons',2)

x_vec = 1:500;

direct_image = Output_counters(:,x_vec,1)/max(Output_counters(:,x_vec,3),[],'all');
direct_and_global_image = Output_counters(:,x_vec,3)/max(Output_counters(:,x_vec,3),[],'all');
global_image = direct_and_global_image-direct_image;

min_val = min([direct_image,direct_and_global_image,global_image],[],'all');
max_val = max([direct_image,direct_and_global_image,global_image],[],'all');

direct_image_norm = imadjust(direct_image,[min_val, max_val],[1,0]);
direct_and_global_image_norm = imadjust(direct_and_global_image,[min_val, max_val],[1,0]);
global_image_norm = direct_image_norm-direct_and_global_image_norm;

write_im_for_pres(direct_image_norm,'Counting Bigger Only Direct Photons',2)
write_im_for_pres(direct_and_global_image_norm,'Counting Bigger Direct and Global Photons',2)
write_im_for_pres(global_image_norm,'Counting Bigger Only Global Photons',2)

%%

x_lims_mat = [ 1, 1508;...
              30,  280;...
               1,  650];

% figure
for ii=1:size(x_lims_mat,1)
    x_vec = x_lims_mat(ii,1):x_lims_mat(ii,2);

    direct_image = Output_counters(:,x_vec,1);
    direct_and_global_image = Output_counters(:,x_vec,3);
    global_image = direct_and_global_image-direct_image;

%     subplot(size(x_lims_mat,1),1,ii)
%     imshow(global_image,[]); colormap('bone');
    
    fprintf('Lim: [%d,%d]. E_global/E_tot = %.2f.\n',x_lims_mat(ii,1),x_lims_mat(ii,2),sum(global_image,'all')/sum(direct_and_global_image,'all'))
%     title(sprintf('E_{global}/E_{tot} = %.2f',sum(global_image,'all')/sum(direct_and_global_image,'all')))
end

%%
slice_number = 1;

figure;
direct_image = Output_counters(:,:,1);
direct_and_global_image = Output_counters(:,:,3);
global_image = direct_and_global_image-direct_image;

shift_val = 222;
shifted_global_image = [global_image(:,(end-shift_val+1):end),global_image(:,1:(end-shift_val))];

limits_image = zeros(size(global_image));
limits_image(:, 30+shift_val) = max(global_image,[],'all');
limits_image(:,280+shift_val) = max(global_image,[],'all');
limits_image(:,754) = max(global_image,[],'all');

imshow(shifted_global_image+limits_image,[]); colormap('bone');
write_im_for_pres(shifted_global_image+limits_image,'Counters All Global',1)

%% Checking spectral scattering:
figure;
for ii=1:6
    subplot(6,1,ii);
    global_image = Output_counters(:,:,5+ii)-Output_counters(:,:,15+ii);
    shift_val = 222;
    shifted_global_image = [global_image(:,(end-shift_val+1):end),global_image(:,1:(end-shift_val))];
    [min_val, max_val] = bounds(shifted_global_image,'all');
    fprintf('Bin %d. Min: %d .Max: %d\n',ii,min_val,max_val);
    imagesc(shifted_global_image); colormap('bone');
    write_im_for_pres(shifted_global_image,sprintf('Global Counters Bin %d Max %d',ii,max_val),1);
end
%%
figure;
global_image_stacked = reshape(permute(Output_counters(:,:,6:11)-Output_counters(:,:,16:21),[2,1,3]),1508,80*6)';
imagesc(global_image_stacked); colormap('bone');