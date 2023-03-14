clear; clc;
used_energy_vec = csvread('210923/used_spectrum.txt');
detected_energy_vec = round(csvread('210923/detected_spectrum.txt'));
%%
[source_energy_vec, ~] = histcounts(used_energy_vec,(0:150));
source_sum_energy = sum(source_energy_vec);
source_energy_vec = source_energy_vec/source_sum_energy;
source_max_energy = max(source_energy_vec);

[detector_energy_vec, ~] = histcounts(detected_energy_vec,(0:150));
detector_sum_energy = sum(detector_energy_vec);
detector_energy_vec = detector_energy_vec/source_sum_energy;
detector_max_energy = max(detector_energy_vec);

%%
num_bins = 6;

% Source equal histogram:
energy_sums = cumsum(source_energy_vec);
energy_per_bin = 1/num_bins;

edges_vec_src_eq_his(1) = 16;
edges_vec_src_eq_his(num_bins+1) = 120;
for ii=2:num_bins
    edges_vec_src_eq_his(ii) = find(energy_sums>(energy_per_bin*(ii-1)),1);
end

edges_vec_src_eq_his = round(edges_vec_src_eq_his);

% Source equal expectation:
energy_sums = cumsum(source_energy_vec.*(1:length(source_energy_vec)));
energy_per_bin = energy_sums(end)/num_bins;

edges_vec_src_eq_exp(1) = 16;
edges_vec_src_eq_exp(num_bins+1) = 120;
for ii=2:num_bins
    edges_vec_src_eq_exp(ii) = find(energy_sums>(energy_per_bin*(ii-1)),1);
end

edges_vec_src_eq_exp = round(edges_vec_src_eq_exp);

% Detector equal histogram:
energy_sums = cumsum(detector_energy_vec);
energy_per_bin = sum(detector_energy_vec)/num_bins;

edges_vec_det_eq_his(1) = 16;
edges_vec_det_eq_his(num_bins+1) = 120;
for ii=2:num_bins
    edges_vec_det_eq_his(ii) = find(energy_sums>(energy_per_bin*(ii-1)),1);
end

edges_vec_det_eq_his = round(edges_vec_det_eq_his);

% Detector equal expectation:
energy_sums = cumsum(detector_energy_vec.*(1:length(detector_energy_vec)));
energy_per_bin = energy_sums(end)/num_bins;

edges_vec_det_eq_exp(1) = 16;
edges_vec_det_eq_exp(num_bins+1) = 120;
for ii=2:num_bins
    edges_vec_det_eq_exp(ii) = find(energy_sums>(energy_per_bin*(ii-1)),1);
end

edges_vec_det_eq_exp = round(edges_vec_det_eq_exp);

edge_vec_mat = [edges_vec_src_eq_his;edges_vec_src_eq_exp;edges_vec_det_eq_his;edges_vec_det_eq_exp];
%%
effective_energy = zeros(size(edge_vec_mat)-[0,1]);
% sums_vec = cumsum(detector_energy_vec.*(1:150));
for ii=1:size(effective_energy,2)
    for jj=1:size(effective_energy,1)
%         effective_energy(jj,ii) = (sums_vec(edge_vec_mat(jj,ii+1))-sums_vec(edge_vec_mat(jj,ii)-1))/(edge_vec_mat(jj,ii+1)-edge_vec_mat(jj,ii));
        effective_energy(jj,ii) = sum(detector_energy_vec(edge_vec_mat(jj,ii):edge_vec_mat(jj,ii+1)).*(edge_vec_mat(jj,ii):edge_vec_mat(jj,ii+1)))/sum(detector_energy_vec(edge_vec_mat(jj,ii):edge_vec_mat(jj,ii+1)));
    end
end
%%

figure;
subplot(3,1,1);
plot(source_energy_vec); hold on;
plot(detector_energy_vec)
title('Energy Distribution')
xlabel('Photon Eneregy [keV]')
ylabel('Probability');
legend('Source','Detector');



for ee=1:size(edge_vec_mat,1)
    
    edges_vec = edge_vec_mat(ee,:)';
    
    effective_energy_locations = 1+edges_vec(1:end-1);
    effective_energy = diff(edges_vec)/2+edges_vec(1:end-1);

    color_mat = [[0 0.4470 0.7410]; ...
                [0.8500 0.3250 0.0980]; ...
                [0.9290 0.6940 0.1250]; ...
                [0.4940 0.1840 0.5560]; ...
                [0.4660 0.6740 0.1880]; ...
                [0.3010 0.7450 0.9330]; ...
                [0.6350 0.0780 0.1840]];

%     figure;
    subplot(3,2,2+ee);
    if ee<3
        energy_vec = source_energy_vec;
        plot(energy_vec); hold on;
        legend_cell{1} = 'Source Energy';
        max_energy = source_max_energy;
    else
        energy_vec = detector_energy_vec;
        plot(energy_vec); hold on;
        legend_cell{1} = 'Detector Energy';
        max_energy = detector_max_energy;
    end
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

%     for ii=1:length(effective_energy_locations)
%         effective_energy(ii) = sum(energy_vec(edges_vec(ii):edges_vec(ii+1)).*((edges_vec(ii):edges_vec(ii+1))))/sum(energy_vec(edges_vec(ii):edges_vec(ii+1)));
%         if mod(ii,2)==1
%             text(effective_energy_locations(ii),0.95*max_energy,sprintf('Effective\n%.2f [keV]',effective_energy(ii)))
%         else
%             text(effective_energy_locations(ii),1.05*max_energy,sprintf('Effective\n%.2f [keV]',effective_energy(ii)))
%         end
%     end
end

%%

save_plot = 0;

h1=figure;
plot(source_energy_vec,'LineWidth',4); hold on;
plot(detector_energy_vec,'LineWidth',4)
if save_plot==0
%     title('Energy Distribution')
    xlabel('Photon Eneregy [keV]')
    ylabel('Intensity');
end
% legend('Source','Detector','','','','','','','','','','','','','','','','','','','','','','','','');
% bin_types_vec = {'Source Histogram - \intI_0dE','Source Expectation - \intEI_0dE','Detector Histogram - \intIdE','Detector Expectation - \intEIdE'};
% bin_types_vec = {[sprintf('Source Histogram\n'),'                      \intI_0dE'],...
%                  [sprintf('Source Expectation\n'),'                      \intEI_0dE'],...
%                  [sprintf('Detector Histogram\n'),'                      \intIdE'],...
%                  [sprintf('Detector Expectation\n'),'                      \intEIdE']};
bin_types_vec = {'Source Histogram','Source Expectation','Detector Histogram','Detector Expectation'};
bin_types_vec = bin_types_vec([1;3;2;4]);

edge_vec_mat = [edges_vec_src_eq_his;edges_vec_src_eq_exp;edges_vec_det_eq_his;edges_vec_det_eq_exp];
edge_vec_mat  = edge_vec_mat([1;3;2;4],:);

color_mat = [[0 0.4470 0.7410]; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]; ...
            [0.3010 0.7450 0.9330]; ...
            [0.6350 0.0780 0.1840]];

% text_x_vec = 91*[1,1,1,1];
text_x_vec = [95,93,93,91];
for ee=1:size(edge_vec_mat,1)
    
    edges_vec = edge_vec_mat(ee,:)';
    
    effective_energy_locations = 1+edges_vec(1:end-1);
    effective_energy = diff(edges_vec)/2+edges_vec(1:end-1);

    for ii=1:(length(edges_vec)-1)
        plot([edges_vec(ii);edges_vec(ii)],[((ee-1)/4)*source_max_energy;(ee/4)*source_max_energy],'Color',color_mat(ii+1,:));
        rectangle('Position', [edges_vec(ii), ((ee-1)/4)*source_max_energy, edges_vec(ii+1)-edges_vec(ii), (1/4)*source_max_energy], ...
                    'FaceColor', [color_mat(ii+1,:), 0.3], ...
                    'EdgeColor', [color_mat(ii+1,:), 0.3]);
        if ii==(length(edges_vec)-1)
            text_x = edge_vec_mat(4,ii);
            text_x = text_x_vec(ee);
            text_y = ((ee-1)/4)*source_max_energy + (1/5)*source_max_energy;
        end
    end
    if save_plot==0
        text(text_x,text_y,bin_types_vec{ee},'FontSize',24)
    end

end
legend('Source','Detector','','','','','','','','','','','','','','','','','','','','','','','','');

set(gca,'FontSize',28);
xticks([16,40:20:120]);
xlim([10, 130]);
% yt = yticks;
yticks([0,0.02,0.04]);

if save_plot==1
    saveas(h1,sprintf('ImagesForPresentations/BinOptions.png'))
    close(h1)
end

%%

save_plot=0;
h1=figure;
plot(source_energy_vec,'LineWidth',4); hold on;
plot(detector_energy_vec,'LineWidth',4)
if save_plot==0
    title('Energy Distribution - 120 [kVp] Source')
    xlabel('Photon Eneregy [keV]')
    ylabel('Probability');
end

color_mat = [[0 0.4470 0.7410]; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]; ...
            [0.3010 0.7450 0.9330]; ...
            [0.6350 0.0780 0.1840]];
        
edges_vec = edge_vec_mat(4,:)';

effective_energy_locations = 1+edges_vec(1:end-1);
effective_energy = diff(edges_vec)/2+edges_vec(1:end-1);

for ii=1:(length(edges_vec)-1)
    plot([edges_vec(ii);edges_vec(ii)],[0;source_max_energy],'Color',color_mat(ii+1,:));
    rectangle('Position', [edges_vec(ii), 0, edges_vec(ii+1)-edges_vec(ii), source_max_energy], ...
                'FaceColor', [color_mat(ii+1,:), 0.2], ...
                'EdgeColor', [color_mat(ii+1,:), 0.2]);
end


legend('Source','Detector','','','','','','');

set(gca,'FontSize',20);
xticks([16,40:20:120]);
xlim([10, 130]);
yt = yticks;
yticks(yt(1:2:end));
if save_plot
    saveas(h1,sprintf('ImagesForPresentations/DetExpBins.png'))
end
