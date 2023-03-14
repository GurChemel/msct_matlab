clear; clc;

file_list = dir('csv_results\*.csv');

for ii=1:length(file_list)
    file_name = file_list(ii).name;
    run_number = sscanf(file_name, 'run_%1d_');
    loss_type  = file_name(7:(end-4));
    Mats.(loss_type)(:,run_number+1) = Import_Loss_From_Csv(['csv_results\',file_name]);
end

figure;
losses = fields(Mats);
for ii=1:length(losses)
    subplot(length(losses),1,ii);
    data_to_plot = (Mats.(losses{ii}));
    plot(data_to_plot)
    title(losses{ii})
    legend({'Multispectral','Normal'})
    ylim([min(0,min(data_to_plot,[],'all')),max(0,max(data_to_plot,[],'all'))])
    xlim([1,size(data_to_plot,1)])
end
