function [linear_atten_3d_mat] = MaterialsAttenuations(enable_plots)

if nargin < 1
    enable_plots = 0;
end

addpath('PhotonAttenuation\');
Materials = load('Materials.mat');
load('Materials.mat','Row_names','Element_Z_Count');
Element_names = Row_names(1:(end-1));

active_elements = {'H','O','P','Ca','I'};
num_active_elements = length(active_elements);
relevant_z_indices = zeros(num_active_elements,1);
for ii=1:num_active_elements
    relevant_z_indices(ii) = find(strcmp(Element_names,active_elements{ii}));
end
relevant_z_count = Element_Z_Count(relevant_z_indices);

% Materials_names = {'Adiposetissue','Water','Iodineblood','dryribwater','Muscle','Blood','Kidney','Liver'};
% Materials_names = {'Adiposetissue','Water','Iodineblood','dryribwater','Muscle','Blood','Liver'};
Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Iodineblood','Redmarrow','Muscle','Liver','dryribwater'};
% Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater'};

num_materials = length(Materials_names);
density_vec = zeros(1,num_materials);
fracs_mat   = zeros(num_active_elements,num_materials);
for ii=1:num_materials
    fracs_mat(:,ii) = Materials.(Materials_names{ii})(relevant_z_indices);
    density_vec(ii) = Materials.(Materials_names{ii})(end);
end
density_mat = repmat(density_vec,num_active_elements,1).*fracs_mat;

%%
% E = (20:2:120)*1e-3;    %in MeV;
% E = [40,60,80,100]*1e-3;    %in MeV;
E_old = [31.31,41.54,50.39,58.94,69.08,90.85]*1e-3;    %in MeV;
E_det_exp = [36.7648   50.4728   59.7230   68.4902   80.8639   98.9679]*1e-3;    %in MeV;
E = E_det_exp([2,4,5,6]);
E_pairs_mat = nchoosek(E,2);
E_triplet = E_det_exp([2,4,6]);

marker_colors_mat = [[0, 0, 1]               ; ...
                     [0, 0.5, 0]             ; ...
                     [1, 0, 0]               ; ...
                     [0, 0.75, 0.75]         ; ...
                     [0.75, 0, 0.75]         ; ...
                     [0.75, 0.75, 0]         ; ...
                     [0, 0.4470, 0.7410]     ; ...
                     [0.8500, 0.3250, 0.0980]; ...
                     [0.9290, 0.6940, 0.1250]; ...
                     [0.4940, 0.1840, 0.5560]; ...
                     [0.4660, 0.6740, 0.1880]; ...
                     [0.3010, 0.7450, 0.9330]; ...
                     [0.6350, 0.0780, 0.1840]];


if enable_plots==1 || enable_plots==5
    figure;
end

linear_atten_3d_mat = zeros(length(Materials_names),2,size(E_pairs_mat,1));

for e_index = 1:size(E_pairs_mat,1)
    mac = PhotonAttenuationQ(relevant_z_count, E_pairs_mat(e_index,:), 'mac')';
    linear_atten = squeeze(sum(repmat(density_mat,1,1,length(E_pairs_mat(e_index,:))).*permute(repmat(mac,1,1,num_materials),[1,3,2]),1));

    linear_atten_3d_mat(:,:,e_index) = linear_atten;
    
    if enable_plots==1
        subplot(2,3,e_index)
        h1=plot(linear_atten(1,1),linear_atten(1,2),'o','Color',marker_colors_mat(1,:),'MarkerFaceColor',marker_colors_mat(1,:),'MarkerSize',15); hold on;
        for ii=2:length(linear_atten(:,1))
            plot(linear_atten(ii,1),linear_atten(ii,2),'o','Color',marker_colors_mat(ii,:),'MarkerFaceColor',marker_colors_mat(ii,:),'MarkerSize',15);
        end
        plot([0,0.6],[0,0.6]);
        xlim([0,0.6])
        ylim([0,0.6])
        title('Mass Attenuation Coefficients in [^{cm^2}/_g] scatter plot')
        xlabel(sprintf('E = %.2f [keV]',E_pairs_mat(e_index,1)*1e3));
        ylabel(sprintf('E = %.2f [keV]',E_pairs_mat(e_index,2)*1e3));
        legend(Materials_names,'Location','southeast')
    end
    if enable_plots==2
        h1=figure;
        plot(linear_atten(1,1),linear_atten(1,2),'o','Color',marker_colors_mat(1,:),'MarkerFaceColor',marker_colors_mat(1,:),'MarkerSize',15); hold on;
        for ii=2:length(linear_atten(:,1))
            plot(linear_atten(ii,1),linear_atten(ii,2),'o','Color',marker_colors_mat(ii,:),'MarkerFaceColor',marker_colors_mat(ii,:),'MarkerSize',15);
        end
        set(gca,'FontSize',15);
        xt = xticks;
        xticks(xt(1:2:end));
        yt = yticks;
        yticks(yt(1:2:end));
        saveas(h1,sprintf('ImagesForPresentations/det_exp_bins/MatAtt_%.2f_%.2f.png',E_pairs_mat(e_index,:)*1e3))
        close(h1)
    end
    if enable_plots==4
        h1=figure;
        plot(linear_atten(1,1),linear_atten(1,2),'o','Color',marker_colors_mat(1,:),'MarkerFaceColor',marker_colors_mat(1,:),'MarkerSize',15); hold on;
        for ii=2:length(linear_atten(:,1))
            plot(linear_atten(ii,1),linear_atten(ii,2),'o','Color',marker_colors_mat(ii,:),'MarkerFaceColor',marker_colors_mat(ii,:),'MarkerSize',15);
        end
        set(gca,'FontSize',15);
        xlim([0,0.4])
        ylim([0,0.4])
        xt = xticks;
        xticks(xt(1:2:end));
        yt = yticks;
        yticks(yt(1:2:end));
        saveas(h1,sprintf('ImagesForPresentations/det_exp_bins_scaled/MatAtt_%.2f_%.2f.png',E_pairs_mat(e_index,:)*1e3))
        close(h1)
        
        figure;
        plot(linear_atten(1,1),linear_atten(1,2),'o','Color',marker_colors_mat(1,:),'MarkerFaceColor',marker_colors_mat(1,:),'MarkerSize',15); hold on;
        for ii=2:length(linear_atten(:,1))
            plot(linear_atten(ii,1),linear_atten(ii,2),'o','Color',marker_colors_mat(ii,:),'MarkerFaceColor',marker_colors_mat(ii,:),'MarkerSize',15);
        end
        title('Mass Attenuation Coefficients in [^{cm^2}/_g] scatter plot')
        xlabel(sprintf('E = %.2f [keV]',E_pairs_mat(e_index,1)*1e3));
        ylabel(sprintf('E = %.2f [keV]',E_pairs_mat(e_index,2)*1e3));
        edited_materials = {'Air', 'Lung Inhale', 'Tissue', 'Breat Mammary', 'Water' ...
                            'Iodine Blode', 'Marrow', 'Muscle', 'Liver', 'Bone'};
        legend(edited_materials,'Location','southeast')
        
    end
    
    if enable_plots==5
        mac = PhotonAttenuationQ(relevant_z_count, [E_det_exp(e_index),0.001:0.001:0.120], 'mac')';
        linear_atten = squeeze(sum(repmat(density_mat,1,1,121).*permute(repmat(mac,1,1,num_materials),[1,3,2]),1));
        
        load('source_energy_vec.mat');
        linear_atten(:,2) = sum(linear_atten(:,2:end).*repmat(source_energy_vec(1:120),num_materials,1),2);
        linear_atten(:,3:end) = [];
        
%         linear_atten_3d_mat(:,:,e_index) = linear_atten;
        subplot(2,3,e_index)
        h1=plot(linear_atten(1,1),linear_atten(1,2),'o','Color',marker_colors_mat(1,:),'MarkerFaceColor',marker_colors_mat(1,:),'MarkerSize',15); hold on;
        for ii=2:length(linear_atten(:,1))
            plot(linear_atten(ii,1),linear_atten(ii,2),'o','Color',marker_colors_mat(ii,:),'MarkerFaceColor',marker_colors_mat(ii,:),'MarkerSize',15);
        end
        plot([0,0.6],[0,0.6]);
        xlim([0,0.6])
        ylim([0,0.6])
        title('Mass Attenuation Coefficients in [^{cm^2}/_g] scatter plot')
        xlabel(sprintf('E = %.2f [keV]',E_det_exp(e_index)*1e3));
%         ylabel(sprintf('E = %.2f [keV]',E_pairs_mat(e_index,2)*1e3));
        ylabel(sprintf('E = Wideband'));
        legend(Materials_names,'Location','southeast')
    end

end


if enable_plots==3
    mac = PhotonAttenuationQ(relevant_z_count, E_triplet, 'mac')';
    linear_atten = squeeze(sum(repmat(density_mat,1,1,length(E_triplet)).*permute(repmat(mac,1,1,num_materials),[1,3,2]),1));

    r = 0.01;
    [x_s,y_s,z_s] = sphere(50);
    figure;
    for ii=1:size(linear_atten,1)
%         x0 = 16.5; y0 = 14.85; z0 = 9.15;
        x = x_s*r + linear_atten(ii,1); %x0;
        y = y_s*r + linear_atten(ii,2); %y0;
        z = z_s*r + linear_atten(ii,3); %z0;
        surf(x,y,z,repmat(permute(marker_colors_mat(ii,:)',[2,3,1]),[size(z),1]),'FaceColor','interp','EdgeColor','none')
        hold on;
    end
end

end