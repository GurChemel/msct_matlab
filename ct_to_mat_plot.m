clear; clc;

save_plot = 0;
LoadGParams;

att_threshold = [gParamStruct.dict_mat_to_dens_values(:,1);2];


color_mat = [[0 0.4470 0.7410]; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]; ...
            [0.3010 0.7450 0.9330]; ...
            0*[0.1000 0.1000 0.1000]; ...
            [0.6350 0.0780 0.1840]];

figure;
for ii=1:(length(att_threshold)-1)
%     [att_threshold(ii),att_threshold(ii+1)-att_threshold(ii)]
    plot([att_threshold(ii)],[0],'Color',color_mat(ii,:),'LineWidth',5); hold on;
    rectangle('Position', [att_threshold(ii), 0, att_threshold(ii+1)-att_threshold(ii), 1], ...
                'FaceColor', [color_mat(ii,:), 1], ...
                'EdgeColor', [color_mat(ii,:), 1]);
end
xlim([-0.1,2.1]);
ylim([0,2]);
if save_plot==0
%     title('Materials Decomposition')
    xlabel('Mass Density [g/cm^3]','FontSize',15)
end

mat_names = {
    'Air'           , ...
    'Lung Inhale'   , ...
    'Adipose Tissue (Fat)', ...
    'Breast Mammary' , ...
    'Water'         , ...
    'Muscle'        , ...
    'Liver'         , ...
    'Bone'};
set(gca,'ytick',[])
legend(mat_names,'Orientation','horizontal','Location','North','FontSize',15);

%%
h1=figure;
plot(source_energy_vec); hold on;
plot(detector_energy_vec)
if save_plot==0
    title('Energy Distribution')
    xlabel('Photon Eneregy [keV]')
    ylabel('Probability');
end
% legend('Source','Detector','','','','','','','','','','','','','','','','','','','','','','','','');
bin_types_vec = {'Source Histogram - \intPdE','Source Expectation - \intPEdE','Detector Histogram - \intPdE','Detector Expectation - \intPEdE'};
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
            text_x = edges_vec(ii+1);
            text_y = ((ee-1)/4)*source_max_energy + (1/8)*source_max_energy;
        end
    end
    if save_plot==0
        text(text_x,text_y,bin_types_vec{ee})
    end

end
legend('Source','Detector','','','','','','','','','','','','','','','','','','','','','','','','');

if save_plot==0
    set(gca,'FontSize',15);
    xticks([16,40:20:120]);
    xlim([10, 130]);
    yt = yticks;
    yticks(yt(1:2:end));
    saveas(h1,sprintf('ImagesForPresentations/BinOptions.png'))
    close(h1)
end