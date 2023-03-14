function [nn] = SmartWeightedNearestNeighbour(data_points,bin_weights,query_points,factor_per_different_label)

%     wighted_data_points = data_points.*repmat(

    query_points_kronned = kron(query_points,ones(size(data_points,1),1));
    data_points_repmated = repmat(data_points,size(query_points,1),1);
    
    bin_weights_repmated = repmat(bin_weights,size(data_points_repmated,1),1);
    
    p = 2;
    diff_per_query = nthroot(sum(bin_weights_repmated.*(abs(query_points_kronned - data_points_repmated).^p),2),p);

    data_diff = reshape(diff_per_query,size(data_points,1),[]);
    [~,first_nn] = min(data_diff,[],1);
    
    first_nn = reshape(first_nn,60,60,[]);
    final_nn = first_nn';
    
    id_relations_norm = [[0.9542, 0.0000, 0.0457, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000];...
                         [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000];...
                         [0.0525, 0.0000, 0.8397, 0.0000, 0.0308, 0.0426, 0.0000, 0.0256, 0.0087, 0.0000];...
                         [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000];...
                         [0.0003, 0.0000, 0.2155, 0.0000, 0.6707, 0.0549, 0.0000, 0.0577, 0.0008, 0.0000];...
                         [0.0000, 0.0000, 0.0903, 0.0000, 0.0166, 0.8672, 0.0000, 0.0034, 0.0224, 0.0000];...
                         [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000];...
                         [0.0000, 0.0000, 0.1226, 0.0000, 0.0396, 0.0078, 0.0000, 0.6837, 0.0050, 0.1413];...
                         [0.0000, 0.0000, 0.1549, 0.0000, 0.0021, 0.1876, 0.0000, 0.0187, 0.6367, 0.0000];...
                         [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0822, 0.0000, 0.9178]];
    id_relations = 1-id_relations_norm;
    id_relations = id_relations-id_relations.*eye(size(id_relations_norm,1));
    
    plot_progress = 0;
    if plot_progress
        figure;
    end
    num_iters = 1;
    for iter_num = 1:num_iters
        label_diff = zeros(size(data_diff));

        for slice=1:size(first_nn,3)
            single_image_material = first_nn(:,:,slice);
            single_image_material_padded = ones(size(single_image_material)+[2,2]);
            single_image_material_padded(2:61,2:61) = single_image_material;
            single_image_material_neighbours(:,1) = reshape(single_image_material_padded(1:60,2:61),[],1);
            single_image_material_neighbours(:,2) = reshape(single_image_material_padded(2:61,1:60),[],1);
            single_image_material_neighbours(:,3) = reshape(single_image_material_padded(3:62,2:61),[],1);
            single_image_material_neighbours(:,4) = reshape(single_image_material_padded(2:61,3:62),[],1);
            single_image_material_neighbours(:,5) = reshape(single_image_material_padded(1:60,1:60),[],1);
            single_image_material_neighbours(:,6) = reshape(single_image_material_padded(3:62,1:60),[],1);
            single_image_material_neighbours(:,7) = reshape(single_image_material_padded(3:62,1:60),[],1);
            single_image_material_neighbours(:,8) = reshape(single_image_material_padded(3:62,3:62),[],1);

%             label_diff_mat = 4*factor_per_different_label*ones(10,60*60);
            label_diff_mat = zeros(size(data_points,1),60*60);
            for ii=1:(60*60)
                for jj=1:size(data_points,1)
                    for kk=1:8
                        label_diff_mat(jj,ii) = label_diff_mat(jj,ii) + factor_per_different_label*id_relations(single_image_material_neighbours(ii,kk),jj);
                    end
%                     label_diff_mat(jj,ii) = label_diff_mat(jj,ii) + factor_per_different_label*id_relations(single_image_material_neighbours(ii,2),jj);
%                     label_diff_mat(jj,ii) = label_diff_mat(jj,ii) + factor_per_different_label*id_relations(single_image_material_neighbours(ii,3),jj);
%                     label_diff_mat(jj,ii) = label_diff_mat(jj,ii) + factor_per_different_label*id_relations(single_image_material_neighbours(ii,4),jj);
                end
            end
            label_diff(:,(1:(60*60))+(slice-1)*(60*60)) = label_diff_mat;
        end

        [~,final_nn] = min(data_diff+label_diff,[],1);
        first_nn = reshape(final_nn,60,60,[]);
        if plot_progress
            subplot(2,ceil(num_iters/2),iter_num);
            imagesc(reshape(final_nn',60,60));
        end
    end
    nn = final_nn';
    
   
end

