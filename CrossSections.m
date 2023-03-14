load('ElementXS.mat')
% %%
% [B,I] = sort(ZVec);
% figure;
% plot(Photon(60,I(1:end-2))); hold on;
% plot(Compton(60,I(1:end-2)));
% plot(Rayleigh(60,I(1:end-2))); hold off;
% figure;
% plot(Photon(20:100,24)); hold on;
% plot(Compton(20:100,24));
% plot(Rayleigh(20:100,24)); hold off;
%%

absorbtion_at_60 = Photon(60,:);
scatter_at_60 = Compton(60,:) + Rayleigh(60,:);

        %     1   2   3   4   5   6   7    8
Names_vec = {'H','C','N','O','P','K','Ca','Cu'};
wanted_z = [1,6,7,8,15,19,20,29];
indices = zeros(1,length(wanted_z));
for ii=1:length(wanted_z)
    indices(ii) = find(wanted_z(ii)==ZVec,1);
end


figure;
subplot(2,2,1)
scatter(scatter_at_60(indices),absorbtion_at_60(indices),'filled');
grid on;
xlabel('\sigma Scatter [cm^2]');
ylabel('\sigma Absorption [cm^2]');
title('Cross Sections at E=60[KeV]')
dx = 1e-25; dy = 1e-24; % displacement so the text does not overlay the data points
text(scatter_at_60(indices)+dx, absorbtion_at_60(indices)+dy, Names_vec, 'Fontsize', 10);

subplot(2,2,2)
scatter(Rayleigh(60,indices),Compton(60,indices),'filled');
grid on;
xlabel('\sigma Rayleigh [cm^2]');
ylabel('\sigma Compton [cm^2]');
title('Cross Sections at E=60[KeV]')
dx = -1e-25; dy = 5e-25; % displacement so the text does not overlay the data points
text(Rayleigh(60,indices)+dx,Compton(60,indices)+dy, Names_vec, 'Fontsize', 10);
%%
figure;
sel_e_vec = [40,60,80,100];
for ii=1:length(sel_e_vec)
    subplot(1,length(sel_e_vec),ii)
    scatter(Rayleigh(sel_e_vec(ii),indices),Compton(sel_e_vec(ii),indices),'filled');
    pbaspect([1,1,1]);
    grid on;
    xlabel('\sigma Rayleigh [cm^2]');
    ylabel('\sigma Compton [cm^2]');
    title(sprintf('Cross Sections at E=%d[KeV]',sel_e_vec(ii)))
    dx = -1e-25; dy = 5e-25; % displacement so the text does not overlay the data points
    text(Rayleigh(sel_e_vec(ii),indices)+dx,Compton(sel_e_vec(ii),indices)+dy, Names_vec, 'Fontsize', 10);
end
%%
figure;
subplot(2,1,1);
semilogy(Rayleigh(15:120,indices));
title('\sigma Rayleigh [cm^2]');
legend(Names_vec);
subplot(2,1,2);
semilogy(Compton(15:120,indices));
title('\sigma Compton [cm^2]');
legend(Names_vec);

%%
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
title('Linear Attenuation');
xlim([0,150])
name_vec{4} = [num2str(alpha),'*',name_vec{1},'+',num2str(beta),'*',name_vec{3}];
legend(name_vec)

% 0.68*O+0.37*Cu = Ca;
o_idx = 4;
ca_idx = 7;
cu_idx = 8;

Na  = 6.02214076e23;
Mu  = 0.999999e-3;
Ao  = 15.999*Mu;
Aca = 40.078*Mu;
Acu = 63.546*Mu;

o_rayleigh  = (Na/Ao)*Rayleigh(15:120,o_idx);
ca_rayleigh = (Na/Aca)*Rayleigh(15:120,ca_idx);
cu_rayleigh = (Na/Acu)*Rayleigh(15:120,cu_idx);
o_compton   = (Na/Ao)*Compton(15:120,o_idx);
ca_compton  = (Na/Aca)*Compton(15:120,ca_idx);
cu_compton  = (Na/Acu)*Compton(15:120,cu_idx);

o_cu_rayleigh = o_rayleigh*0.68+cu_rayleigh*0.37;
o_cu_compton  = o_compton*0.68+cu_compton*0.37;

% figure;
% subplot(1,3,2);
nexttile; 
semilogy([o_rayleigh,ca_rayleigh,cu_rayleigh,o_cu_rayleigh],'LineWidth',3); hold on;
ylabel('\sigma [cm^2]');
xlabel('Energy [KeV]');
legend({'O','Ca','Cu','0.68*O+0.37*Cu'});
title('Rayleigh Cross-Section');
% subplot(1,3,3);
nexttile; 
semilogy([o_compton,ca_compton,cu_compton,o_cu_compton],'LineWidth',3);
ylabel('\sigma [cm^2]');
xlabel('Energy [KeV]');
legend({'O','Ca','Cu','0.68*O+0.37*Cu'});
title('Compton Cross-Section');

%%

rayleigh_at_60 = Rayleigh(60,indices);

angle_vec = linspace(-pi,pi,200);
phase_vector = 0.375*(1+cos(angle_vec).^2);%.*sin(angle_vec);

rayleigh_per_phase = phase_vector'*rayleigh_at_60;
figure;
subplot(2,1,1);
plot(phase_vector)
subplot(2,1,2);
plot(repmat(angle_vec',1,length(indices)),rayleigh_per_phase)
legend(Names_vec)