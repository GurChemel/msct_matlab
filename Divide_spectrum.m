clear; clc;
energy_vec = csvread('spectrum.txt');
max_energy = max(energy_vec);


num_bins = 6;
divide_type = 'equal_histogram';
% divide_type = 'equal_expectation';

if strcmp(divide_type,'equal_histogram')
    energy_sums = cumsum(energy_vec);
    energy_per_bin = 1/num_bins;

    edges_vec(1) = 16;
    edges_vec(num_bins+1) = 120;
    for ii=2:num_bins
        edges_vec(ii) = find(energy_sums>(energy_per_bin*(ii-1)),1);
    end

    edges_vec = round(edges_vec);
end
if strcmp(divide_type,'equal_expectation')
    energy_sums = cumsum(energy_vec'.*(1:length(energy_vec)));
    energy_per_bin = energy_sums(end)/num_bins;

    edges_vec(1) = 16;
    edges_vec(num_bins+1) = 120;
    for ii=2:num_bins
        edges_vec(ii) = find(energy_sums>(energy_per_bin*(ii-1)),1);
    end

    edges_vec = round(edges_vec);
end

effective_energy_locations = 1+edges_vec(1:end-1);
effective_energy = diff(edges_vec)/2+edges_vec(1:end-1);

figure;
plot(energy_vec);
title('Source Output Energy Distribution - 120 [kVp]')
xlabel('Photon Eneregy [keV]')
ylabel('Probability');

color_mat = [[0 0.4470 0.7410]; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]; ...
            [0.3010 0.7450 0.9330]; ...
            [0.6350 0.0780 0.1840]];

figure;
plot(energy_vec); hold on;
legend_cell{1} = 'Source Energy';
for ii=1:(length(edges_vec)-1)
    plot([edges_vec(ii);edges_vec(ii)],[0;max_energy],'Color',color_mat(ii+1,:));
    legend_cell{ii+1} = sprintf('Bin %d',ii);
    rectangle('Position', [edges_vec(ii), 0, edges_vec(ii+1)-edges_vec(ii), max_energy], ...
                'FaceColor', [color_mat(ii+1,:), 0.3], ...
                'EdgeColor', [color_mat(ii+1,:), 0.3]);
end
legend(legend_cell)
title('Source Output Energy Distribution - 120 [kVp]')
xlabel('Photon Eneregy [keV]')
ylabel('Probability');

for ii=1:length(effective_energy_locations)
    effective_energy(ii) = sum(energy_vec(edges_vec(ii):edges_vec(ii+1)).*(edges_vec(ii):edges_vec(ii+1))')/sum(energy_vec(edges_vec(ii):edges_vec(ii+1)));
    if mod(ii,2)==1
        text(effective_energy_locations(ii),0.95*max_energy,sprintf('Effective\n%.2f [keV]',effective_energy(ii)))
    else
        text(effective_energy_locations(ii),1.05*max_energy,sprintf('Effective\n%.2f [keV]',effective_energy(ii)))
    end
end
