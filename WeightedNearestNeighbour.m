function [nn] = WeightedNearestNeighbour(data_points,bin_weights,query_points)

%     wighted_data_points = data_points.*repmat(

    query_points_kronned = kron(query_points,ones(size(data_points,1),1));
    data_points_repmated = repmat(data_points,size(query_points,1),1);
    
    bin_weights_repmated = repmat(bin_weights,size(data_points_repmated,1),1);
    
    diff_per_query = sqrt(sum(bin_weights_repmated.*((query_points_kronned - data_points_repmated).^2),2));

    [~,nn] = min(reshape(diff_per_query,size(data_points,1),[]),[],1);
    
    nn = nn';
end

