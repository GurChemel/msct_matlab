clear; clc;

solid_gt_1 = load('220222/GT_solid_1.mat');
solid_gt_1 = solid_gt_1.I_GT;
% solid_vol_1 = load('220215/xcat_reduced_solid_1.mat');
% solid_vol_1_id = solid_vol_1.xcat_id;
% solid_vol_1_density = solid_vol_1.xcat_density;

solid_gt_2 = load('220222/GT_solid_2.mat');
solid_gt_2 = solid_gt_2.I_GT;
% solid_vol_2 = load('220215/xcat_reduced_solid_2.mat');
% solid_vol_2_id = solid_vol_2.xcat_id;
% solid_vol_2_density = solid_vol_2.xcat_density;

tmp_solid_gt_1 = cat(2,solid_gt_1(:,1600:end,:),solid_gt_1(:,1:200,:));
tmp_solid_gt_2 = cat(2,solid_gt_2(:,1600:end,:),solid_gt_2(:,1:200,:));
figure
for ii=1:11
    subplot(11,2,2*ii-1);
    imagesc(squeeze(tmp_solid_gt_1(ii,:,:))');% colormap('bone');
    if ii==1
        title('[Solid 1]');
    end
    subplot(11,2,2*ii);
    imagesc(squeeze(tmp_solid_gt_2(ii,:,:))');% colormap('bone');
    if ii==1
        title('[Solid 2]');
    end
end

tmp_solid_gt_1 = cat(2,solid_gt_1(:,1770:end,20:60),solid_gt_1(:,1:30,20:60));
tmp_solid_gt_2 = cat(2,solid_gt_2(:,1770:end,20:60),solid_gt_2(:,1:30,20:60));
figure
for ii=1:11
%     subplot(11,2,2*ii-1);
%     imagesc(squeeze(solid_gt_1(ii,:,:))');% colormap('bone');
%     subplot(11,2,2*ii);
%     imagesc(squeeze(solid_gt_2(ii,:,:))');% colormap('bone');
    subplot(11,2,2*ii-1);
    imagesc([squeeze(tmp_solid_gt_1(ii,:,:))',squeeze(tmp_solid_gt_2(ii,:,:))']);% colormap('bone');
    if ii==1
        title('[Solid 1, Solid 2]');
    end
    subplot(11,3,3*ii);
    imagesc([squeeze(tmp_solid_gt_1(ii,:,:))'-squeeze(tmp_solid_gt_2(ii,:,:))']);% colormap('bone');
    if ii==1
        title('[Solid 1 - Solid 2]');
    end
end
%%
E = linspace(10e-3,150e-3,400);

% z_vec = [6,8,20];
% name_vec = {'C','O','Ca'};

z_vec = [8,20,29];
name_vec = {'O','Ca','Cu'};

mac = PhotonAttenuationQ(z_vec, E, 'mac');

% dens_vec = [15.999, 40.078, 63.54];
dens_vec = [1, 1, 1];
mac_linear = mac .* repmat(dens_vec,size(mac,1),1);

alpha_beta = round((mac_linear([20,190],[1,3])^-1)*mac_linear([20,190],2),2);
alpha = 0.68;
beta = 0.37;
total = 1;
% alpha = alpha_beta(1);
% beta = alpha_beta(2);

mac_linear(:,4) = total*(mac_linear(:,1)*alpha+mac_linear(:,3)*beta);
figure;
% subplot(1,3,1);
tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
nexttile; 
semilogy(repmat(E'*1e3,1,size(mac_linear,2)),(mac_linear),'LineWidth',3)
xlabel('Energy [KeV]');
ylabel('\mu [1/cm]');
xlim([0,150])
name_vec{4} = [num2str(alpha),'*',name_vec{1},'+',num2str(beta),'*',name_vec{3}];
legend(name_vec)

%%
E = linspace(10e-3,150e-3,400);

z_vec = [1,5,6,7,8,20,26];
% name_vec = {'O','Ca','Cu','I'};
% z_vec = [20,53,64,83];
% name_vec = {'O','I','Gb','Bi'};

mac = PhotonAttenuationQ(z_vec, E, 'mac');

figure
semilogy(repmat(E'*1e3,1,size(mac,2)),(mac))
xlabel('E [kev]')
xlim([0,150])
ylim([0.1 1000])
legend()
set(gca,'FontSize',24);