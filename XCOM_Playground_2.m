clear;
clc;

files = dir('XCOM_Data/*.csv');
for ii=1:length(files)
    atom_number = files(ii).name(1:(end-4));
    if str2double(atom_number) < 21 || str2double(atom_number) == 53
        eval(['atom_',atom_number,'=csvread(''XCOM_Data/',atom_number,'.csv'',2,0);']);
        eval(['atom_',atom_number,'(:,1)=atom_',atom_number,'(:,1)*1000;']);
    end
end
atom_20 = atom_20([1:5,8:19],:);
atom_53 = atom_53([1,5:8,14,17:22,25:29],1:6);

               %           H,     C,     N,     O,     P,   Ca,     I 
blood_frac_density   = [0.102,  0.11, 0.033, 0.745, 0.001,    0,     0]*1.06;
tissue_frac_density  = [0.114, 0.598, 0.007, 0.278,     0,    0,     0]*0.95;
bone_frac_density    = [0.034, 0.155, 0.042, 0.435, 0.103,0.225,     0]*1.55;
blood_iodine_density = [0.100, 0.109, 0.033, 0.739, 0.001,    0, 0.008]*1.09;

atom_data = zeros(17,6,6);
atom_data(:,:,1) = atom_1;
atom_data(:,:,2) = atom_6;
atom_data(:,:,3) = atom_7;
atom_data(:,:,4) = atom_8;
atom_data(:,:,5) = atom_14;
atom_data(:,:,6) = atom_20;
atom_data(:,:,7) = atom_53;

materials_data = zeros(17,6,4);
materials_data(:,:,1) = cat(2,atom_data(:,1,1),sum(atom_data(:,2:end,:).*repmat(permute(tissue_frac_density,[1,3,2]),17,5,1),3));
materials_data(:,:,2) = cat(2,atom_data(:,1,1),sum(atom_data(:,2:end,:).*repmat(permute(blood_frac_density,[1,3,2]),17,5,1),3));
materials_data(:,:,3) = cat(2,atom_data(:,1,1),sum(atom_data(:,2:end,:).*repmat(permute(bone_frac_density,[1,3,2]),17,5,1),3));
materials_data(:,1,4) = materials_data(:,1,1);
materials_data(:,:,5) = cat(2,atom_data(:,1,1),sum(atom_data(:,2:end,:).*repmat(permute(blood_iodine_density,[1,3,2]),17,5,1),3));
name_vec = {'Soft Tissue \rho=0.95[g/cm^3]','Blood \rho=1.06[g/cm^3]','Bone \rho=1.55[g/cm^3]'};

%%
           %     Rayleigh              Compton
legend_vec = {'Coherent Scatter.','Incoher. Scatter.','Photoel. Absorb.','Tot. w/ Coherent','Tot. wo/ Coherent'};

rayleigh_index = 2;
compton_index = 3;
pa_index = 4;
total_index = 5;

figure;
for ii=1:3
    subplot(1,3,ii);
    semilogy(repmat(materials_data(:,1,ii),1,5),materials_data(:,2:6,ii));
    legend(legend_vec)
    title(name_vec{ii});
end

%%


           %     Rayleigh              Compton
%legend_vec = {'Coherent Scatter.','Incoher. Scatter.','Photoel. Absorb.','Tot. w/ Coherent','Tot. wo/ Coherent'};
legend_vec = {'Coherent Scatter.','Incoher. Scatter.','Photoel. Absorb.','Tot. w/ Coherent'};

rayleigh_index = 2;
compton_index = 3;
pa_index = 4;
total_index = 5;

axis_font_size=20;
for ii=[1,3]
    figure('DefaultAxesFontSize',axis_font_size)
%     subplot(1,3,ii);
%     set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',axis_font_size);
    semilogy(repmat(materials_data(:,1,ii),1,4),materials_data(:,2:5,ii),'LineWidth',3);
    ylim([1e-3,1e4]);
%     legend(legend_vec)
%     title(name_vec{ii});
%     set(gca,'YTickLabel',get(gca,'YTickLabel'),'fontsize',axis_font_size);
end
%%
figure;
ind_to_plot=[5];
ii=1;
semilogy(repmat(materials_data(:,1,ii),1,length(ind_to_plot)),materials_data(:,ind_to_plot,ii),'LineWidth',3); hold on;
ii=3;
semilogy(repmat(materials_data(:,1,ii),1,length(ind_to_plot)),materials_data(:,ind_to_plot,ii),'--','LineWidth',3); hold on;
ylim([1e-3,1e4]);


%%
alpha = 1.04;
beta = 0.045;

materials_data(:,2:6,4) = materials_data(:,2:6,1)*alpha+materials_data(:,2:6,3)*beta;
name_vec{4} = ['Soft Tissue \rho=',num2str(alpha*0.95),'[g/cm^3] + Bone \rho=',num2str(round(beta*1.55,2)),'[g/cm^3]'];

ylim_vec = [1e-3,1e1];

figure;
tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
nexttile; 
semilogy(materials_data(:,1,1),materials_data(:,pa_index,1),'LineWidth',3); hold on;
semilogy(materials_data(:,1,2),materials_data(:,pa_index,2),'LineWidth',3);
semilogy(materials_data(:,1,3),materials_data(:,pa_index,3),'LineWidth',3);
semilogy(materials_data(:,1,4),materials_data(:,pa_index,4),'LineWidth',3);
xlabel({'Energy [KeV]';'(a)'});
ylabel('\sigma [1/cm]');
title('Photon Absorbtion Cross Section');
xlim([15,100])
% ylim(ylim_vec)
legend(name_vec)

nexttile; 
semilogy(materials_data(:,1,1),materials_data(:,rayleigh_index,1),'LineWidth',3); hold on;
semilogy(materials_data(:,1,2),materials_data(:,rayleigh_index,2),'LineWidth',3);
semilogy(materials_data(:,1,3),materials_data(:,rayleigh_index,3),'LineWidth',3);
semilogy(materials_data(:,1,4),materials_data(:,rayleigh_index,4),'LineWidth',3);
ylabel('\sigma [1/cm]');
xlabel({'Energy [KeV]';'(b)'});
xlim([15,100])
% ylim(ylim_vec)
legend(name_vec);
title('Rayleigh Cross-Section');

% nexttile; 
% semilogy(materials_data(:,1,1),materials_data(:,compton_index,1),'LineWidth',3); hold on;
% semilogy(materials_data(:,1,2),materials_data(:,compton_index,2),'LineWidth',3);
% semilogy(materials_data(:,1,3),materials_data(:,compton_index,3),'LineWidth',3);
% semilogy(materials_data(:,1,4),materials_data(:,compton_index,4),'LineWidth',3);
% ylabel('\sigma [1/cm]');
% xlabel({'Energy [KeV]';'(c)'});
% xlim([15,100])
% % ylim(ylim_vec)
% legend(name_vec);
% title('Compton Cross-Section');

nexttile; 
semilogy(materials_data(:,1,1),materials_data(:,total_index,1),'LineWidth',3); hold on;
semilogy(materials_data(:,1,2),materials_data(:,total_index,2),'LineWidth',3);
semilogy(materials_data(:,1,3),materials_data(:,total_index,3),'LineWidth',3);
semilogy(materials_data(:,1,4),materials_data(:,total_index,4),'LineWidth',3);
ylabel('\sigma [1/cm]');
xlabel({'Energy [KeV]';'(d)'});
xlim([15,100])
% ylim(ylim_vec)
legend(name_vec);
title('All Cross-Sections');

%%

% name_vec = {'Soft Tissue \rho=0.95[g/cm^3]','Blood \rho=1.06[g/cm^3]','Bone \rho=1.55[g/cm^3]'};
% name_vec{4} = ['Soft Tissue \rho=',num2str(alpha*0.95),'[g/cm^3] + Bone \rho=',num2str(round(beta*1.55,2)),'[g/cm^3]'];
name_vec = {'Soft Tissue','Blood','Bone','Soft Tissue + Bone'};

font_size = 20;
axis_font_size = 22;
figure;
% tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
semilogy(materials_data(:,1,1),materials_data(:,pa_index,1),'-','LineWidth',3); hold on;
semilogy(materials_data(:,1,2),materials_data(:,pa_index,2),':','LineWidth',3);
semilogy(materials_data(:,1,3),materials_data(:,pa_index,3),'-.','LineWidth',3);
semilogy(materials_data(:,1,4),materials_data(:,pa_index,4),'--','LineWidth',3);
% title('Photon Absorbtion Cross Section','FontSize',font_size);
xlim([15,100])
yticks([0.01 1 100])
yticklabels({'0.01','1','100'})
set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',axis_font_size);
xlabel('Energy [KeV]','FontSize',font_size);
ylabel('Photon absorption coefficient [1/cm]','FontSize',font_size);
%set(gca,'YTickLabel',get(gca,'YTickLabel'),'fontsize',font_size);

% ylim(ylim_vec)
legend(name_vec,'FontSize',font_size)

figure; 
semilogy(materials_data(:,1,1),materials_data(:,rayleigh_index,1),'-','LineWidth',3); hold on;
semilogy(materials_data(:,1,2),materials_data(:,rayleigh_index,2),':','LineWidth',3);
semilogy(materials_data(:,1,3),materials_data(:,rayleigh_index,3),'-.','LineWidth',3);
semilogy(materials_data(:,1,4),materials_data(:,rayleigh_index,4),'--','LineWidth',3);
xlim([15,100])
yticks([0.01 0.1 1])
yticklabels({'0.01','0.1','1'})
% ylim(ylim_vec)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',axis_font_size);
set(gca,'YTickLabel',get(gca,'YTickLabel'),'fontsize',axis_font_size);
ylabel('Rayleigh scattering coefficient [1/cm]','FontSize',font_size);
xlabel('Energy [KeV]','FontSize',font_size);
legend(name_vec,'FontSize',font_size);
% title('Rayleigh Cross-Section','FontSize',font_size);

figure; 
semilogy(materials_data(:,1,1),materials_data(:,compton_index,1),'-','LineWidth',3); hold on;
semilogy(materials_data(:,1,2),materials_data(:,compton_index,2),':','LineWidth',3);
semilogy(materials_data(:,1,3),materials_data(:,compton_index,3),'-.','LineWidth',3);
semilogy(materials_data(:,1,4),materials_data(:,compton_index,4),'--','LineWidth',3);
xlim([15,100])
% ylim(ylim_vec)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',axis_font_size);
set(gca,'YTickLabel',get(gca,'YTickLabel'),'fontsize',axis_font_size);
ylabel('Compton scattering coefficient [1/cm]','FontSize',font_size);
xlabel('Energy [KeV]','FontSize',font_size);
legend(name_vec,'FontSize',font_size);
% title('Rayleigh Cross-Section','FontSize',font_size);

figure; 
semilogy(materials_data(:,1,1),materials_data(:,total_index,1),'-','LineWidth',3); hold on;
semilogy(materials_data(:,1,2),materials_data(:,total_index,2),':','LineWidth',3);
semilogy(materials_data(:,1,3),materials_data(:,total_index,3),'-.','LineWidth',3);
semilogy(materials_data(:,1,4),materials_data(:,total_index,4),'--','LineWidth',3);
xlim([15,100])
yticks([0.1 1 10])
yticklabels({'0.1','1','10'})
set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',axis_font_size);
set(gca,'YTickLabel',get(gca,'YTickLabel'),'fontsize',axis_font_size);
% ylim(ylim_vec)
ylabel('Attenuation coefficient [1/cm]','FontSize',font_size);
xlabel('Energy [KeV]','FontSize',font_size);
legend(name_vec,'FontSize',font_size);
% title('All Cross-Sections','FontSize',font_size);


%%
color_mat = [[0 0.4470 0.7410]; ...
            [0.8500 0.3250 0.0980]; ...
            [0.9290 0.6940 0.1250]; ...
            [0.4940 0.1840 0.5560]; ...
            [0.4660 0.6740 0.1880]; ...
            [0.3010 0.7450 0.9330]; ...
            [0.6350 0.0780 0.1840]];

% color_mat = [[0,0,1];...
%              [0,1,0];...
%              [0,1,1];...
%              [1,0,0];...
%              [1,0,1];...
%              [1,1,0];...
%              [0,0,0];...
%              [1,1,1]]*.8;
figure;
% semilogy(materials_data(:,1,2),materials_data(:,pa_index,2),'-.','LineWidth',3,'Color',color_mat(1,:)); hold on;
% semilogy(materials_data(:,1,4),materials_data(:,pa_index,4),'-.','LineWidth',3,'Color',color_mat(2,:));
% semilogy(materials_data(:,1,2),materials_data(:,rayleigh_index,2),'--','LineWidth',3,'Color',color_mat(1,:));
% semilogy(materials_data(:,1,4),materials_data(:,rayleigh_index,4),'--','LineWidth',3,'Color',color_mat(2,:));
% semilogy(materials_data(:,1,2),materials_data(:,compton_index,2),':','LineWidth',3,'Color',color_mat(1,:));
% semilogy(materials_data(:,1,4),materials_data(:,compton_index,4),':','LineWidth',3,'Color',color_mat(2,:));
% semilogy(materials_data(:,1,2),materials_data(:,total_index,2),'LineWidth',3,'Color',color_mat(1,:));
% semilogy(materials_data(:,1,4),materials_data(:,total_index,4),'LineWidth',3,'Color',color_mat(2,:));

semilogy(materials_data(:,1,2),materials_data(:,pa_index,2),'LineWidth',3,'Color',color_mat(1,:)); hold on;
semilogy(materials_data(:,1,4),materials_data(:,pa_index,4),'LineWidth',3,'Color',color_mat(2,:));
semilogy(materials_data(:,1,2),materials_data(:,rayleigh_index,2),'LineWidth',3,'Color',color_mat(3,:));
semilogy(materials_data(:,1,4),materials_data(:,rayleigh_index,4),'LineWidth',3,'Color',color_mat(4,:));
% semilogy(materials_data(:,1,2),materials_data(:,compton_index,2),'LineWidth',3,'Color',color_mat(5,:));
% semilogy(materials_data(:,1,4),materials_data(:,compton_index,4),'LineWidth',3,'Color',color_mat(6,:));
semilogy(materials_data(:,1,2),materials_data(:,total_index,2),'LineWidth',3,'Color',color_mat(5,:));
semilogy(materials_data(:,1,4),materials_data(:,total_index,4),'LineWidth',3,'Color',color_mat(6,:));


xlabel('Energy [KeV]');
ylabel('\sigma [1/cm]');
% title('Photon Absorbtion Cross Section');
% title('Rayleigh Cross-Section');
% title('Compton Cross-Section');
% title('All Cross-Sections');
xlim([15,100])
% ylim(ylim_vec)
% legend(name_vec)

%%

rayleigh_index = 2;
compton_index = 3;
pa_index = 4;
total_index = 5;

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

% materials_data

linear_atten = zeros(3,3,4);
for figure_index=2:5
    cnt=1;
    for ii=[2,3,5]
        linear_atten(cnt,1,figure_index-1) = materials_data(14,figure_index,ii);
        linear_atten(cnt,2,figure_index-1) = materials_data(16,figure_index,ii);
        linear_atten(cnt,3,figure_index-1) = materials_data(17,figure_index,ii);
        cnt = cnt + 1;
    end
end
%%
max_val_vec = [0.1,0.3,0.4,0.65];
titles_vec = {'Rayleigh Scattering','Compton Scattering','Photoelectric Absorption','Total'};

% figure;
for jj=1:4
%     switch jj
%         case 1
%             subplot(1,2,1);
%             shape_dot = 'o';
%         case 2
%             shape_dot = 'd';
%         case 3
%             subplot(1,2,2);
%             shape_dot = 'o';
%         case 4
%             shape_dot = 'd';
%     end
    h1=figure;
    yyaxis left
    plot(linear_atten(1,3,jj),linear_atten(1,1,jj),shape_dot,'Color',marker_colors_mat(1,:),'MarkerFaceColor',marker_colors_mat(1,:),'MarkerSize',15); hold on;
    for ii=2:length(linear_atten(:,1))
        plot(linear_atten(ii,3,jj),linear_atten(ii,1,jj),shape_dot,'Color',marker_colors_mat(ii,:),'MarkerFaceColor',marker_colors_mat(ii,:),'MarkerSize',15);
    end
    pbaspect([1,1,1]);
    max_val = max_val_vec(jj);
    xlim([0,max_val])
    ylim([0,max_val])
    xlabel('E = 100 [keV]');
    ylabel('E = 50 [keV]');

    yyaxis right
    plot(linear_atten(1,3,jj),linear_atten(1,2,jj),shape_dot,'Color',marker_colors_mat(3+1,:),'MarkerFaceColor',marker_colors_mat(3+1,:),'MarkerSize',15); hold on;
    for ii=2:length(linear_atten(:,1))
        plot(linear_atten(ii,3,jj),linear_atten(ii,2,jj),shape_dot,'Color',marker_colors_mat(3+ii,:),'MarkerFaceColor',marker_colors_mat(3+ii,:),'MarkerSize',15);
    end
    xlim([0,max_val])
    ylim([0,max_val])
    ylabel('E = 80 [keV]');
    plot([0,max_val],[0,max_val],'k-');

    legend({'Blood Left Y-Axis','Bone Left Y-Axis','Blood-Iodine Left Y-Axis',...
            'Blood Right Y-Axis','Bone Right Y-Axis','Blood-Iodine Right Y-Axis'},'Location','southeast');
    title(titles_vec{jj});
    file_name_str = titles_vec{jj};
    file_name_str(file_name_str==' ')='_';
    
    saveas(h1,sprintf('ImagesForPresentations/cross_sec_dots/cross_section_%s.png',file_name_str))
    close(h1)
end
%%
figure;
e_vec_for_now = [50,80,100];
for e_off=0:1
    subplot(1,2,1+e_off)
    shape_vec = {'+','x','*','o'};
    for fig_idx=4
        plot(linear_atten(1,1+e_off,fig_idx),linear_atten(1,2+e_off,fig_idx),shape_vec{fig_idx},'Color',marker_colors_mat(1,:),'MarkerFaceColor',marker_colors_mat(1,:),'MarkerSize',15); hold on;
        for ii=2:length(linear_atten(:,1))
            plot(linear_atten(ii,1+e_off,fig_idx),linear_atten(ii,2+e_off,fig_idx),shape_vec{fig_idx},'Color',marker_colors_mat(ii,:),'MarkerFaceColor',marker_colors_mat(ii,:),'MarkerSize',15);
        end
    end
    max_val = 0.65;
    plot([0,max_val],[0,max_val]);
    xlim([0,max_val])
    ylim([0,max_val])
%     legend({'blood rayleigh','bone rayleigh','blood iodine rayleigh',...
%             'blood compton','bone compton','blood iodine compton',...
%             'blood absorption','bone absorption','blood iodine absorption',...
%             'blood total','bone total','blood iodine total'})
    legend({'Blood','Bone','Blood-Iodine'})
    pbaspect([1,1,1]);
    set(gca,'FontSize',15);
    title('Attenuation Coefficients in [^{1}/_{cm}] scatter plot')
    xlabel(sprintf('E = %.2f [keV]',e_vec_for_now(e_off+1)));
    ylabel(sprintf('E = %.2f [keV]',e_vec_for_now(e_off+2)));
    %         ylabel(sprintf('E = %.2f [keV]',E_pairs_mat(e_index,2)*1e3));
    % ylabel(sprintf('E = Wideband'));
end