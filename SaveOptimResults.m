clear; clc;

folder_name = '220325';
files = dir([folder_name,'/*.csv']);

num_step = size(readmatrix([folder_name,'/',files(1).name]),1);
num_files = 12;

if length(files)~=num_files
    error('Number of files error');
end

OptimResults = zeros(num_step,num_files);
for file_id=1:num_files
    file_data = readmatrix([folder_name,'/',files(file_id).name]);
%     disp(['Working on: ',files(file_id).name]);
    run_type = [];
    if contains(files(file_id).name,'accum_scat')
        run_type = 4;
    elseif contains(files(file_id).name,'accum_lin')
        run_type = 3;
    elseif contains(files(file_id).name,'per_bin_scat')
        run_type = 2;
    elseif contains(files(file_id).name,'per_bin_lin')
        run_type = 1;
    end
    loss_type = [];
    if contains(files(file_id).name,'loss')
        loss_type = 0;
    elseif contains(files(file_id).name,'epsilon')
        loss_type = 1;
    elseif contains(files(file_id).name,'delta')
        loss_type = 2;
    end
    OptimResults(:,loss_type*4+run_type) = file_data(:,3);
end

save(['OptimResults_',folder_name,'.mat'],'OptimResults');
