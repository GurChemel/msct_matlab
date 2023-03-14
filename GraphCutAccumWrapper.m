function [labels] = GraphCutAccumWrapper(input_image,smothness_factor)

    if nargin < 2
        smothness_factor = 0.03;
    end

    Materials_names = {'Air','Water','Muscle','dryribwater','Redmarrow','Iodineblood','Adiposetissue'};
    energy_kev = 60;
    [linear_atten_mat, ~] = PerEnergyMaterialsAttenuations(Materials_names, energy_kev);

    query_points = reshape(input_image,[],size(input_image,3));

    query_points_kronned = kron(query_points,ones(size(linear_atten_mat,1),1));
    data_points_repmated = repmat(linear_atten_mat,size(query_points,1),1);

    p = 2;
    diff_per_query = nthroot((abs(query_points_kronned - data_points_repmated).^p),p);

    DataCost = reshape(reshape(diff_per_query,size(linear_atten_mat,1),[])',size(input_image,1),size(input_image,2),[]);

    id_relations_norm = [[0.9542, 0.0457, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000];...
                         [0.0525, 0.8397, 0.0308, 0.0426, 0.0256, 0.0087, 0.0000];...
                         [0.0003, 0.2155, 0.6707, 0.0549, 0.0577, 0.0008, 0.0000];...
                         [0.0000, 0.0903, 0.0166, 0.8672, 0.0034, 0.0224, 0.0000];...
                         [0.0000, 0.1226, 0.0396, 0.0078, 0.6837, 0.0050, 0.1413];...
                         [0.0000, 0.1549, 0.0021, 0.1876, 0.0187, 0.6367, 0.0000];...
                         [0.0000, 0.0000, 0.0000, 0.0000, 0.0822, 0.0000, 0.9178]];

    SmoothnessCost = 1-(id_relations_norm+id_relations_norm')/2;
    SmoothnessCost = SmoothnessCost-SmoothnessCost.*eye(7);

    [gch] = GraphCut('open', DataCost, smothness_factor*SmoothnessCost);

    [gch, labels] = GraphCut('swap', gch);

    [gch] = GraphCut('close', gch);
    
    labels = labels+1; % Starting from 1 for use as indices.
end