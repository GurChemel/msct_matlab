clear;clc;

files = dir('outCount\run_0outputCounter*.csv');

for ii=1:length(files)
    bin_idx = str2double(files(ii).name(20:(end-4)))-14;
    outputCnt(:,:,bin_idx) = import_runOutputs(['outCount\',files(ii).name]);
end

files = dir('outCount\run_0outputEnergyDep*.csv');

for ii=1:length(files)
    bin_idx = str2double(files(ii).name(22:(end-4)))-14;
    outputEnergy(:,:,bin_idx) = import_runOutputs(['outCount\',files(ii).name]);
end
%%

energy_centers = [31.3090   41.5364   50.3864   58.9392   69.0782   90.8517].';

figure;
for selected_bin = 1:6
    subplot(4,6,selected_bin);
    imagesc(outputCnt(10:310,:,selected_bin).');
    title('Real Counter'); colorbar()
    subplot(4,6,selected_bin+6);
    imagesc(outputEnergy(10:310,:,selected_bin).');
    title('Real Energy'); colorbar()
    subplot(4,6,selected_bin+12);
    imagesc(energy_centers(selected_bin)*outputCnt(10:310,:,selected_bin).');
    title('Estimated Energy'); colorbar()
    subplot(4,6,selected_bin+18);
    max_energy = max(abs(outputEnergy(10:310,:,selected_bin).'),[],'all');
    max_diff = max(abs(outputEnergy(10:310,:,selected_bin).' - energy_centers(selected_bin)*outputCnt(10:310,:,selected_bin).'),[],'all');
    imagesc(outputEnergy(10:310,:,selected_bin).' - energy_centers(selected_bin)*outputCnt(10:310,:,selected_bin).');
    title(sprintf('Max Diff: %.2f %%',max_diff/max_energy)); colorbar()
end
%%
N_vec = squeeze(sum(sum(outputCnt,1),2))';
fprintf('N_vec = [%.0f, %.0f, %.0f, %.0f, %.0f, %.0f]\n',N_vec)
Sqrt_N_vec = sqrt(N_vec);

estimated_sqrt_n_vec = [0.9291, 0.9685, 0.9791, 0.9837, 0.9855, 0.9872];

mat = [Sqrt_N_vec/Sqrt_N_vec(6);estimated_sqrt_n_vec/estimated_sqrt_n_vec(6)]
