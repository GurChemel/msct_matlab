function [] = FindSolids()

E = linspace(10e-3,150e-3,400);

z_vec = [8,20,29];
name_vec = {'O','Ca','Cu'};

mac = PhotonAttenuationQ(z_vec, E, 'mac');

alpha_init = 0.68;
beta_init = 0.37;
total_init = 1;
% 
% mac(:,4) = total*(mac(:,1)*alpha+mac(:,3)*beta);
% 
% semilogy(repmat(E'*1e3,1,size(mac,2)),(mac),'.')
% xlabel('E [kev]')
% xlim([0,150])
% name_vec{4} = [num2str(alpha),'*',name_vec{1},'+',num2str(beta),'*',name_vec{3}];
% legend(name_vec)

%prepare figure and guidata struct
h=struct;
h.f=figure;
h.ax=axes('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.1 0.1 0.6 0.8]);
h.slider=uicontrol('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.7 0.1 0.1 0.8],...
    'Style','Slider',...
    'BackgroundColor',[1 1 1],...
    'Min',0,'Max',1,'Value',alpha_init,'SliderStep',[0.01,0.1],...
    'Callback',@sliderCallback);
h.slider2=uicontrol('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.8 0.1 0.1 0.8],...
    'Style','Slider',...
    'BackgroundColor',[1 1 1],...
    'Min',0,'Max',1,'Value',beta_init,'SliderStep',[0.01,0.1],...
    'Callback',@sliderCallback2);
h.slider3=uicontrol('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.9 0.1 0.1 0.8],...
    'Style','Slider',...
    'BackgroundColor',[1 1 1],...
    'Min',0,'Max',5,'Value',total_init,'SliderStep',[0.01,0.1],...
    'Callback',@sliderCallback3);
%store image database to the guidata struct as well
h.database=mac;
h.alpha = alpha_init;
h.beta = beta_init;
h.total = total_init;
guidata(h.f,h)
%trigger a callback
sliderCallback(h.slider)
function sliderCallback(hObject,~)
    
    h.alpha = get(hObject,'Value');
    alpha = h.alpha;
    beta = h.beta;
    total = h.total;
    mac=guidata(hObject).database;

    mac(:,4) = total*(mac(:,1)*alpha+mac(:,3)*beta);
    semilogy(repmat(E'*1e3,1,size(mac,2)),(mac),'.')
    title(sprintf('\\alpha = %.2f. \\beta = %.2f. Total = %.2f',alpha,beta,total));
    xlabel('E [kev]')
    xlim([0,150])
    name_vec{4} = [num2str(alpha),'*',name_vec{1},'+',num2str(beta),'*',name_vec{3}];
    legend(name_vec)
end
function sliderCallback2(hObject,~)
    
    h.beta = get(hObject,'Value');
    alpha = h.alpha;
    beta = h.beta;
    total = h.total;
    mac=guidata(hObject).database;

    mac(:,4) = total*(mac(:,1)*alpha+mac(:,3)*beta);
    semilogy(repmat(E'*1e3,1,size(mac,2)),(mac),'.')
    title(sprintf('\\alpha = %.2f. \\beta = %.2f. Total = %.2f',alpha,beta,total));
    xlabel('E [kev]')
    xlim([0,150])
    name_vec{4} = [num2str(alpha),'*',name_vec{1},'+',num2str(beta),'*',name_vec{3}];
    legend(name_vec)
end
function sliderCallback3(hObject,~)
    
    h.total = get(hObject,'Value');
    alpha = h.alpha;
    beta = h.beta;
    total = h.total;
    mac=guidata(hObject).database;

    mac(:,4) = total*(mac(:,1)*alpha+mac(:,3)*beta);
    semilogy(repmat(E'*1e3,1,size(mac,2)),(mac),'.')
    title(sprintf('\\alpha = %.2f. \\beta = %.2f. Total = %.2f',alpha,beta,total));
    xlabel('E [kev]')
    xlim([0,150])
    name_vec{4} = [num2str(alpha),'*',name_vec{1},'+',num2str(beta),'*',name_vec{3}];
    legend(name_vec)
end
end