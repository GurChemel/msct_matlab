% function [A] = BuildSparseA()

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
date_str = '220519'; %datestr(now,'YYmmDD')
for sourceIndex = 1:Nsources
    curr_source_id  = src_id(sourceIndex)+1;
    curr_source_loc = src_loc(curr_source_id,:);
    if isfile(sprintf('sparseA_%s/source_id_%03d.mat',date_str,curr_source_id))
        fprintf('Skipping Source Number %3d/%d. Source ID %3d. .mat file exists.\n',sourceIndex, Nsources,curr_source_id);
        continue;
    end
    tic
    fprintf('Source Number %3d/%d. Source ID %3d. Location: (%8.3f,%8.3f,%8.3f). Working   0%%',sourceIndex, Nsources, curr_source_id, curr_source_loc);
    SparseA = sparse(Tpixels , Tvoxels);

%     L = [voxelSizeX*nx, voxelSizeY*ny, voxelSizeZ*nz] ;  % your cube dimensions 
%     O = -L/2 ;       % Get the origin of cube so that P is at center 
%     plotcube(L,O,.2,[1 0 0]);   % use function plotcube 
%     hold on

    % for all pixels in detector
    for Row = 1: det_height
        for Column = 1: det_width
            fprintf('\b\b\b\b%3d%%',round(100*((Row-1)*det_width+Column)/(det_width*det_height)));
            pixelIndex = sub2ind( [ det_height , det_width ] , Row , Column );
            
            curr_det_id = det_id_per_source(sourceIndex,pixelIndex);
            mod_curr_det_id = mod(curr_det_id,80);
            if (mod_curr_det_id==79)
                vert_next_det = curr_det_id-1;
            else
                vert_next_det = curr_det_id+1;
            end
            if curr_det_id>=143920
                horz_next_det = curr_det_id-80;
            else
                horz_next_det = curr_det_id+80;
            end
            
            curr_det_loc = det_loc(1+curr_det_id,:);
            next_det_vert_loc = det_loc(1+vert_next_det,:);
            next_det_horz_loc = det_loc(1+horz_next_det,:);
            curr_det_loc_rand = curr_det_loc + (rand-0.5)*(next_det_vert_loc-curr_det_loc) + (rand-0.5)*(next_det_horz_loc-curr_det_loc);
            
            ray_direction = curr_det_loc_rand-curr_source_loc;
            
            theta = acos(ray_direction(3)/norm(ray_direction));
            phi   = acos(ray_direction(1)/norm(ray_direction(1:2)));
            
            dr = 0.1;
            dx = sin(theta)*cos(phi);
            dy = sin(theta)*sin(phi);
            dz = cos(theta);
            
            direction = [dx, dy, dz];
            step_direction = dr*direction;

            [flag ,tmin] = rayBoxIntersection(curr_source_loc, ray_direction, boundry_min, boundry_max);
            if flag == 0
                continue
            end
            
            % Find first interaction with the volume
%             step_mult = 1;
%             max_steps = norm(ray_direction)/dr;
%             while (any((curr_source_loc+step_mult*step_direction) < boundry_min) || any((curr_source_loc+step_mult*step_direction) > boundry_max))
%                 step_mult = step_mult + 1;
%                 if step_mult>max_steps
%                     break
%                 end
%             end
%             if step_mult>max_steps
%                 continue
%             end
% 
%             curr_intersect_loc = (curr_source_loc+step_mult*step_direction);
%             curr_intersect_indices = floor(curr_intersect_loc./voxel_size_vec+(volume_res_vec/2))+1;
%             curr_tmax_vec = curr_intersect_loc-((curr_intersect_indices-1-(volume_res_vec/2)).*voxel_size_vec);
%             if (any(curr_tmax_vec>voxel_size_vec))
%                 error('tMax Calculation bigger than voxel size!')
%             end

            if (tmin<0)
                tmin = 0;
            end
            start   = curr_source_loc + tmin*ray_direction;
            boxSize = boundry_max-boundry_min;

            
            rounded_eps = round((start-boundry_min)./boxSize,5);
            
            x = floor( rounded_eps(1)*nx )+1;
            y = floor( rounded_eps(2)*ny )+1;
            z = floor( rounded_eps(3)*nz )+1;               
            if (x==(nx+1));  x=x-1;  end
            if (y==(ny+1));  y=y-1;  end          
            if (z==(nz+1));  z=z-1;  end

           % calculate step orintation
           if (ray_direction(1)>=0)
                tVoxelX = (x)/nx;
                stepX = 1;
            else
                tVoxelX = (x-1)/nx;
                stepX = -1;  
           end

            if (ray_direction(2)>=0)
                tVoxelY = (y)/ny;
                stepY = 1;
            else
                tVoxelY = (y-1)/ny;
                stepY = -1;
            end

            if (ray_direction(3)>=0)
                tVoxelZ = (z)/nz; 
                stepZ = 1;
            else
                tVoxelZ = (z-1)/nz;
                stepZ = -1;  
            end

            voxelMaxX  = boundry_min(1) + tVoxelX*boxSize(1);
            voxelMaxY  = boundry_min(2) + tVoxelY*boxSize(2);
            voxelMaxZ  = boundry_min(3) + tVoxelZ*boxSize(3);
            tMaxX      = tmin + (voxelMaxX-start(1))/direction(1);
            tMaxY      = tmin + (voxelMaxY-start(2))/direction(2);
            tMaxZ      = tmin + (voxelMaxZ-start(3))/direction(3);
        
            tDeltaX    = voxelSizeX/abs(direction(1));
            tDeltaY    = voxelSizeY/abs(direction(2));
            tDeltaZ    = voxelSizeZ/abs(direction(3));
            
%             tMaxX = curr_tmax_vec(1);
%             tMaxY = curr_tmax_vec(2);
%             tMaxZ = (voxelSizeZ/voxelSizeY)*curr_tmax_vec(3);

            index_x = x;
            index_y = y;
            index_z = z;
            

            ref = tmin; %%0 ;
            Vx = zeros(1,Tvoxels);
            Vy = zeros(1,Tvoxels);
            Vz = zeros(1,Tvoxels);
            RayInVoxel_lengths = zeros(1,Tvoxels);
            curr_v_ray_vec_idx = 1;

            Vx(curr_v_ray_vec_idx) = index_x;
            Vy(curr_v_ray_vec_idx) = index_y;
            Vz(curr_v_ray_vec_idx) = index_z;

            while ( true )

                % add new value
                VoxInd = sub2ind([nx , ny , nz],index_x,index_y,index_z);

                if (tMaxX < tMaxY)
                    if (tMaxX < tMaxZ)
                        RayInVoxel_lengths(curr_v_ray_vec_idx) = tMaxX - ref  ;
                        index_x = index_x + stepX;
                        if (index_x > nx || index_x<1)
                            break;
                        end
                        ref =  tMaxX;
                        tMaxX = tMaxX + tDeltaX;

                    else
                        RayInVoxel_lengths(curr_v_ray_vec_idx) = tMaxZ - ref ;
                        index_z = index_z + stepZ;
                        if (index_z > nz || index_z<1)
                            break;
                        end
                        ref = tMaxZ ;
                        tMaxZ = tMaxZ + tDeltaZ;

                    end
                else
                    if (tMaxY < tMaxZ)
                        RayInVoxel_lengths(curr_v_ray_vec_idx) = tMaxY - ref ;
                        index_y = index_y + stepY;
                        if (index_y > ny || index_y<1)
                            break;
                        end
                        ref =  tMaxY ;
                        tMaxY = tMaxY + tDeltaY;

                    else
                        RayInVoxel_lengths(curr_v_ray_vec_idx) = tMaxZ - ref ;
                        index_z = index_z + stepZ;
                        if (index_z > nz || index_z<1)
                            break;
                        end
                        ref =  tMaxZ ;
                        tMaxZ = tMaxZ + tDeltaZ;

                    end
                end

                if (index_y < 1)
                    disp('ERROR Occurred in CalculateRayLengthInVoxel( )');
                    break ;
                end

                curr_v_ray_vec_idx = curr_v_ray_vec_idx + 1;
                Vx(curr_v_ray_vec_idx) = index_x;
                Vy(curr_v_ray_vec_idx) = index_y;
                Vz(curr_v_ray_vec_idx) = index_z;


            end

            Vx((curr_v_ray_vec_idx+1):end) = [];
            Vy((curr_v_ray_vec_idx+1):end) = [];
            Vz((curr_v_ray_vec_idx+1):end) = [];
            RayInVoxel_lengths((curr_v_ray_vec_idx+1):end) = [];
            
            FullResVoxIndex = [Vx ; Vy ; Vz ];% was befor
            FullResVoxIndex = sub2ind([ nx , ny , nz] , FullResVoxIndex(1,:) , FullResVoxIndex(2,:) , FullResVoxIndex(3,:) );
            
            SparseA(pixelIndex,FullResVoxIndex) = RayInVoxel_lengths ;
%             SparseA = sparse(pixelIndex, FullResVoxIndex, RayInVoxel_lengths, Tpixels, Tvoxels);

%             mat = [curr_det_loc;curr_source_loc];
%             plot3(mat(:,1),mat(:,2),mat(:,3))
%             plot3(first_intersect_loc(1),first_intersect_loc(2),first_intersect_loc(3),'x')
                
            
        end
    end
    
    save(sprintf('sparseA_%s/source_id_%03d.mat',date_str,curr_source_id),'SparseA')
    fprintf('. '); toc
end

%%

% P = [0,0,0] ;   % you center point 
% L = [voxelSizeX*nx, voxelSizeY*ny, voxelSizeZ*nz] ;  % your cube dimensions 
% O = -L/2 ;       % Get the origin of cube so that P is at center 
% plotcube(L,O,.2,[1 0 0]);   % use function plotcube 
% hold on
% 
% mat = [curr_det_loc;curr_source_loc];
% plot3(mat(:,1),mat(:,2),mat(:,3))
% end