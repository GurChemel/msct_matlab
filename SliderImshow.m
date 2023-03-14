function [] = SliderImshow(images)

%prepare figure and guidata struct
h=struct;
h.f=figure;
h.ax=axes('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.1 0.1 0.6 0.8]);
h.slider=uicontrol('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.8 0.1 0.1 0.8],...
    'Style','Slider',...
    'BackgroundColor',[1 1 1],...
    'Min',1,'Max',size(images,3),'Value',1,...
    'Callback',@sliderCallback);
%store image database to the guidata struct as well
h.database=images;
guidata(h.f,h)
%trigger a callback
sliderCallback(h.slider)
function sliderCallback(hObject,~)
    h=guidata(hObject);
    slide_num=round(get(hObject,'Value'));
    imshow(h.database(:,:,slide_num),[]);
    title(['Slide Number: ',num2str(slide_num)]);
end
end