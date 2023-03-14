function [] = FindSolids(images)

alpha_init = min(images,[],'all');
beta_init = max(images,[],'all');
total_init = 1;

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
    'Min',alpha_init,'Max',beta_init,'Value',alpha_init,'SliderStep',[0.01,0.1],...
    'Callback',@sliderCallback);
h.slider2=uicontrol('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.8 0.1 0.1 0.8],...
    'Style','Slider',...
    'BackgroundColor',[1 1 1],...
    'Min',alpha_init,'Max',beta_init,'Value',beta_init,'SliderStep',[0.01,0.1],...
    'Callback',@sliderCallback2);
h.slider3=uicontrol('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.9 0.1 0.1 0.8],...
    'Style','Slider',...
    'BackgroundColor',[1 1 1],...
    'Min',1,'Max',size(images,3),'Value',total_init,...'SliderStep',[0.01,0.1],...
    'Callback',@sliderCallback3);
%store image database to the guidata struct as well
h.database = images;
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
    images_tmp=guidata(hObject).database;
    imshow(images_tmp(:,:,total),[alpha,beta]);
    title(sprintf('Min = %.2f. Max = %.2f. Slide = %d from %d',alpha,beta,total,size(images_tmp,3)));
end
function sliderCallback2(hObject,~)
    
    h.beta = get(hObject,'Value');
    alpha = h.alpha;
    beta = h.beta;
    total = h.total;
    images_tmp=guidata(hObject).database;
    imshow(images_tmp(:,:,total),[alpha,beta]);
    title(sprintf('Min = %.2f. Max = %.2f. Slide = %d from %d',alpha,beta,total,size(images_tmp,3)));
end
function sliderCallback3(hObject,~)
    
    h.total = round(get(hObject,'Value'));
    alpha = h.alpha;
    beta = h.beta;
    total = h.total;
    images_tmp=guidata(hObject).database;
    imshow(images_tmp(:,:,total),[alpha,beta]);
    title(sprintf('Min = %.2f. Max = %.2f. Slide = %d from %d',alpha,beta,total,size(images_tmp,3)));
end
end