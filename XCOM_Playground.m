clear;
clc;

files = dir('XCOM_Data/*.csv');
for ii=1:length(files)
    atom_number = files(ii).name(1:(end-4));
    eval(['atom_',atom_number,'=csvread(''XCOM_Data/',atom_number,'.csv'',2,0);']);
    eval(['atom_',atom_number,'(:,1)=atom_',atom_number,'(:,1)*1000;']);
end

atom_29 = atom_29([1,5:11,14:22],:);
atom_20 = atom_20([1:5,8:19],:);
atom_53 = atom_53([1,5:8,14,17:22,25:29],:);
atom_26 = atom_26([1:7,10:19],:);


atom_6 = atom_26;
atom_7 = atom_20;

% vec = [1:7,10:19];
% sprintf('%5d ',1000*atom_20(:,1))
% sprintf('%5d ',1000*atom_26(:,1))
% sprintf('%5d ',vec)

% %%
%            %     Rayleigh              Compton
% legend_vec = {'Coherent Scatter.','Incoher. Scatter.','Photoel. Absorb.','Tot. w/ Coherent','Tot. wo/ Coherent'};
% 
% 
% figure;
% subplot(1,3,1);
% semilogy(repmat(atom_6(:,1),1,5),atom_6(:,2:6));
% legend(legend_vec)
% title('8 - Oxygen (O)');
% subplot(1,3,2);
% semilogy(repmat(atom_7(:,1),1,5),atom_7(:,2:6))
% legend(legend_vec)
% title('20 - Calcium (Ca)');
% subplot(1,3,3);
% semilogy(repmat(atom_8(:,1),1,5),atom_8(:,2:6))
% legend(legend_vec)
% title('26 - Iron (Fe)');

%%
% In Blood:
O_density = 1.06*0.745;
Fe_density = 1.06*0.001;


name_vec = {'N','C','O'};

alpha = 0.5;
beta = 0.5;
total = 1;



figure;
tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
nexttile; 
semilogy(atom_6(:,1),atom_6(:,4),'LineWidth',3); hold on;
semilogy(atom_7(:,1),atom_7(:,4),'LineWidth',3);
semilogy(atom_8(:,1),atom_8(:,4),'LineWidth',3);
semilogy(atom_8(:,1),alpha*atom_6(:,4)+beta*atom_8(:,4),'LineWidth',3);
xlabel('Energy [KeV]');
ylabel('\mu [1/cm]');
% title('Linear Attenuation');
title('Photon Absorbtion');
xlim([15,100])
name_vec{4} = [num2str(alpha),'*',name_vec{1},'+',num2str(beta),'*',name_vec{3}];
legend(name_vec)

% 0.68*O+0.37*Fe = Ca;
n_idx = 4;
c_idx = 7;
o_idx = 8;

Na  = 1;%6.02214076e23;
Mu  = 1;%0.999999e-3;
An  = 1;%15.999*Mu;
Ac = 1;%40.078*Mu;
Ao = 1;%63.546*Mu;

n_rayleigh  = (Na/An)*atom_6(:,2);
c_rayleigh = (Na/Ac)*atom_7(:,2);
o_rayleigh = (Na/Ao)*atom_8(:,2);
n_compton   = (Na/An)*atom_6(:,3);
c_compton  = (Na/Ac)*atom_7(:,3);
o_compton  = (Na/Ao)*atom_8(:,3);

n_o_rayleigh = n_rayleigh*alpha+o_rayleigh*beta;
n_o_compton  = n_compton*alpha+o_compton*beta;


nexttile; 
semilogy(atom_6(:,1),n_rayleigh,'LineWidth',3); hold on;
semilogy(atom_6(:,1),c_rayleigh,'LineWidth',3);
semilogy(atom_6(:,1),o_rayleigh,'LineWidth',3);
semilogy(atom_6(:,1),n_o_rayleigh,'LineWidth',3);
ylabel('\sigma [cm^2]');
xlabel('Energy [KeV]');
xlim([15,100])
legend(name_vec);
title('Rayleigh Cross-Section');

nexttile; 
semilogy(atom_6(:,1),n_compton,'LineWidth',3); hold on;
semilogy(atom_6(:,1),c_compton,'LineWidth',3);
semilogy(atom_6(:,1),o_compton,'LineWidth',3);
semilogy(atom_6(:,1),n_o_compton,'LineWidth',3);
ylabel('\sigma [cm^2]');
xlabel('Energy [KeV]');
xlim([15,100])
legend(name_vec);
title('Compton Cross-Section');


nexttile; 
semilogy(atom_6(:,1),atom_6(:,5)+n_rayleigh+n_compton,'LineWidth',3); hold on;
semilogy(atom_6(:,1),atom_7(:,5)+c_rayleigh+c_compton,'LineWidth',3);
semilogy(atom_6(:,1),atom_8(:,5)+o_rayleigh+o_compton,'LineWidth',3);
semilogy(atom_6(:,1),alpha*atom_6(:,5)+beta*atom_8(:,5)+n_o_rayleigh+n_o_compton,'LineWidth',3);
ylabel('\sigma [cm^2]');
xlabel('Energy [KeV]');
xlim([15,100])
legend(name_vec);
title('All Cross-Sections');

%%

% rayleigh_sum = [sum(o_rayleigh),sum(c_rayleigh),sum(o_rayleigh),sum(o_o_rayleigh)];
% compton_sum  = [sum(o_compton),sum(c_compton),sum(o_compton),sum(o_o_compton)];
% 
% plot(rayleigh_sum,compton_sum,'.')

% %%
% color_a = [     0, 0.4470, 0.7410];
% color_b = [0.8500, 0.3250, 0.0980];
% color_c = [0.9290, 0.6940, 0.1250];
% color_d = [0.4940, 0.1840, 0.5560];
% 
% line_width = 1.5;
% 
% figure;
% semilogy(atom_6(:,1),atom_6(:,4),'-.','LineWidth',line_width,'Color',color_a); hold on;
% semilogy(atom_20(:,1),atom_20(:,4),'-.','LineWidth',line_width,'Color',color_b);
% semilogy(atom_29(:,1),atom_29(:,4),'-.','LineWidth',line_width,'Color',color_c);
% semilogy(atom_29(:,1),alpha*atom_6(:,5)+beta*atom_29(:,5),'-.','LineWidth',3,'Color',color_d);
% xlabel('Energy [KeV]');
% ylabel('\sigma [cm^2]');
% % title('Linear Attenuation');
% title('Photon Absorbtion');
% xlim([15,100])
% name_vec{4} = [num2str(alpha),'*',name_vec{1},'+',num2str(beta),'*',name_vec{3}];
% legend(name_vec)
% 
% semilogy(atom_6(:,1),o_rayleigh,'--','LineWidth',line_width,'Color',color_a);
% semilogy(atom_6(:,1),c_rayleigh,'--','LineWidth',line_width,'Color',color_b);
% semilogy(atom_6(:,1),cu_rayleigh,'--','LineWidth',line_width,'Color',color_c);
% semilogy(atom_6(:,1),o_cu_rayleigh,'--','LineWidth',line_width,'Color',color_d);
% 
% semilogy(atom_6(:,1),o_compton,':','LineWidth',line_width,'Color',color_a);
% semilogy(atom_6(:,1),c_compton,':','LineWidth',line_width,'Color',color_b);
% semilogy(atom_6(:,1),cu_compton,':','LineWidth',line_width,'Color',color_c);
% semilogy(atom_6(:,1),o_cu_compton,':','LineWidth',line_width,'Color',color_d);
% 
% semilogy(atom_6(:,1),atom_6(:,5)+o_rayleigh+o_compton,'LineWidth',3,'Color',color_a);
% semilogy(atom_6(:,1),atom_20(:,5)+c_rayleigh+c_compton,'LineWidth',3,'Color',color_b);
% semilogy(atom_6(:,1),atom_29(:,5)+cu_rayleigh+cu_compton,'LineWidth',3,'Color',color_c);
% semilogy(atom_6(:,1),alpha*atom_6(:,5)+beta*atom_29(:,5)+o_cu_rayleigh+o_cu_compton,'LineWidth',3,'Color',color_d);
