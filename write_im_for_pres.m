function [] = write_im_for_pres(image,title_str,zoom_factor)
    
    if nargin < 3
        zoom_factor = 4;
    end

    image = image/max(image,[],'all');
    image(image<0) = 0;

    file_name_str = title_str;
    file_name_str(file_name_str==' ')='_';
    file_name_str = ['ImagesForPresentations/',file_name_str,'.png'];
    imwrite(kron(image,ones(zoom_factor)),file_name_str);
end

