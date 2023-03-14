clear; clc;


src_loc = csvread('source_and_det_2\src_loc.csv');
src_orient = csvread('source_and_det_2\src_orient.csv');
det_loc = csvread('source_and_det_2\det_loc.csv');
cut_frames = csvread('source_and_det_2\cut_frames_2.csv');

src_id = cut_frames(1:180,1);
det_id_borders = cut_frames(1:180,[2,3])+1;
det_angle_per_source = zeros(length(src_id),401);
for loop_num=1:length(src_id)
    if diff(det_id_borders(loop_num,:))>0
        det_angle_per_source(loop_num,:) = det_id_borders(loop_num,1):det_id_borders(loop_num,2);
    else
        det_angle_per_source(loop_num,:) = [det_id_borders(loop_num,1):1800,1:det_id_borders(loop_num,2)];
    end
end

det_angle_per_source = det_angle_per_source-1;
det_id_per_source = kron(det_angle_per_source,ones(1,80))*80+repmat(0:79',180,401);

%%
Nsources = size(det_id_per_source,1);

det_height = 80;
det_width  = 401;
Tpixels = size(det_id_per_source,2);

% medium resolution (voxels)
nz = 14;  % vertical_cells;
nx = 256; % horizontal_cells;
ny = 256; % horizontal_cells;
volume_res_vec = [nx, ny, nz];
Tvoxels = nx*ny*nz;

voxelSizeX = 2*0.24;
voxelSizeY = 2*0.24;
voxelSizeZ = 2*1;
voxel_size_vec = [voxelSizeX,voxelSizeY,voxelSizeZ];

boundry_min = -(volume_res_vec.*voxel_size_vec)/2;
boundry_max = +(volume_res_vec.*voxel_size_vec)/2;

%%
num_sources = 2;
figure;
for sourceIndex = 1:num_sources
    curr_source_id  = src_id(sourceIndex)+1;
    curr_source_loc = src_loc(curr_source_id,:);

    origin = curr_source_loc;

    load(sprintf('sparseA_220503/source_id_%03d.mat',curr_source_id));
    [pixel_idx_vec,voxel_idx_vec,ray_length_vec] = find(SparseA);
    subplot(1,num_sources,sourceIndex)
    
    %%
    det_id = 1+(cut_frames(sourceIndex,[2,3])*80);

    if diff(det_id)>0
        plot_det_id = det_id(1):det_id(2);
    else
        plot_det_id = [det_id(1):1800,1:det_id(2)];
    end
    
    scatter3(src_loc(:,1),src_loc(:,2),src_loc(:,3),'.'); hold on;
%     scatter3(det_loc(plot_det_id,1),det_loc(plot_det_id,2),det_loc(plot_det_id,3),'.');

    %%
    det_to_plot = 5;
    selected_pixels = randsample(unique(pixel_idx_vec),det_to_plot)';
    for pixelIndex=selected_pixels(1:det_to_plot)

        relevant_indices = find(pixel_idx_vec==pixelIndex);
        [volX,volY,volZ] = ind2sub([nx,ny,nz],voxel_idx_vec(relevant_indices));
        relevant_ray_length_vec = ray_length_vec(relevant_indices);

        curr_det_loc = det_loc(1+det_id_per_source(sourceIndex,pixelIndex),:);
        direction = curr_det_loc-curr_source_loc;

        hold on;
        text(origin(1), origin(2), origin(3), 'origin');
        plot3(origin(1), origin(2), origin(3), 'k.', 'MarkerSize', 15);
        quiver3(origin(1), origin(2), origin(3), direction(1), direction(2), direction(3));

        vmin = boundry_min;
        vmax = boundry_max;
        BoxVertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
        BoxFaces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
        h = patch('Vertices',BoxVertices,'Faces',BoxFaces,'FaceColor','yellow');
        set(h, 'FaceAlpha', 0.1);
        view(60,30);
    %     axis tight;
        xlim([-280,280])
        ylim([-280,280])
        zlim([ -40, 40])
        xlabel('x');
        ylabel('y');
        zlabel('z');
        grid on;


        for corr_voxel_id = 1:length(relevant_ray_length_vec)
            vmin = (boundry_min + ([volX(corr_voxel_id),volY(corr_voxel_id),volZ(corr_voxel_id)]-1).*voxel_size_vec);
            vmax = (boundry_min + ([volX(corr_voxel_id),volY(corr_voxel_id),volZ(corr_voxel_id)]).*voxel_size_vec);
            smallBoxVertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
            smallBoxFaces    = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];

            h = patch('Vertices', smallBoxVertices, 'Faces', smallBoxFaces, 'FaceColor', 'blue', 'EdgeColor', 'blue');
            set(h,'FaceAlpha',ray_length_vec(corr_voxel_id));
        end
    end
end