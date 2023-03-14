function [labels] = GraphCut3dPerBinWrapper(input_image, smothness_factor, dect)

    if nargin < 2
        smothness_factor = 0.03;
    end
    if nargin < 3
        dect = 0;
    end

    Materials_names = {'Air','Water','Muscle','dryribwater','Redmarrow','Iodineblood','Adiposetissue'};
    
    if dect==1
        active_materials = [1, 4, 7];
    elseif dect==2
        active_materials = [1, 3, 4, 7];
    elseif dect==3
        active_materials = [1, 2, 3, 4, 5];
    else
        active_materials = 1:7;
    end
    Materials_names = Materials_names(active_materials);
    [linear_atten_mat, ~] = PerBinMaterialsAttenuations(Materials_names);
%     bin_weights = [0, 0.2, 0.4, 0.6, 0.8, 1];
    bin_weights = [0, 1, 0.5, 0, 0, 0];

    query_points = reshape(input_image,[],6);

    query_points_kronned = kron(query_points,ones(size(linear_atten_mat,1),1));
    data_points_repmated = repmat(linear_atten_mat,size(query_points,1),1);

    bin_weights_repmated = repmat(bin_weights,size(data_points_repmated,1),1);

    p = 2;
    diff_per_query = nthroot(sum(bin_weights_repmated.*(abs(query_points_kronned - data_points_repmated).^p),2),p);

    DataCost = reshape(reshape(diff_per_query,size(linear_atten_mat,1),[])',size(input_image,1),size(input_image,2),size(input_image,3),[]);

    id_relations_norm = [[0.9542, 0.0457, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000];...
                         [0.0525, 0.8397, 0.0308, 0.0426, 0.0256, 0.0087, 0.0000];...
                         [0.0003, 0.2155, 0.6707, 0.0549, 0.0577, 0.0008, 0.0000];...
                         [0.0000, 0.0903, 0.0166, 0.8672, 0.0034, 0.0224, 0.0000];...
                         [0.0000, 0.1226, 0.0396, 0.0078, 0.6837, 0.0050, 0.1413];...
                         [0.0000, 0.1549, 0.0021, 0.1876, 0.0187, 0.6367, 0.0000];...
                         [0.0000, 0.0000, 0.0000, 0.0000, 0.0822, 0.0000, 0.9178]];

    id_relations_norm = id_relations_norm(active_materials,active_materials);
                     
    SmoothnessCost = 1-(id_relations_norm+id_relations_norm')/2;
    SmoothnessCost = SmoothnessCost-SmoothnessCost.*eye(length(Materials_names));

    [gch] = GraphCut('open', DataCost, smothness_factor*SmoothnessCost);

    [gch, labels] = GraphCut('swap', gch);

    [gch] = GraphCut('close', gch);
    
    labels = labels+1; % Starting from 1 for use as indices.
    
    labels = reshape(labels,size(input_image,1),size(input_image,2),size(input_image,3));
end