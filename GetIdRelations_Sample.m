clear;clc;
load('211005/xcat_reduced.mat','xcat_id');
xcat_materials = load('xcat_materials.mat');

Materials_names = {'Air','Adipose_tissue','Water','Muscle','dry_rib_water','Blood','Red_marrow'};
num_mats = length(Materials_names);

relevant_ids = zeros(1,num_mats);
xcat_shorten_ids = xcat_id;
for ii = 1:num_mats
    relevant_ids(ii) = xcat_materials.(Materials_names{ii});
    xcat_shorten_ids(xcat_id==relevant_ids(ii))=ii;
end
%%
id_relations = zeros(num_mats);
for row=2:59
    for col=2:59
        for slice=1:10
            cur_id = xcat_shorten_ids(row,col,slice);
            if row > 1
                id_relations(cur_id,xcat_shorten_ids(row-1,col,slice)) = id_relations(cur_id,xcat_shorten_ids(row-1,col,slice)) + 1;
            end
            if row < 60
                id_relations(cur_id,xcat_shorten_ids(row+1,col,slice)) = id_relations(cur_id,xcat_shorten_ids(row+1,col,slice)) + 1;
            end
            if col > 1
                id_relations(cur_id,xcat_shorten_ids(row,col-1,slice)) = id_relations(cur_id,xcat_shorten_ids(row,col-1,slice)) + 1;
            end
            if col < 60
                id_relations(cur_id,xcat_shorten_ids(row,col+1,slice)) = id_relations(cur_id,xcat_shorten_ids(row,col+1,slice)) + 1;
            end
        end
    end
end

%%
id_relations_norm = id_relations./sum(id_relations);
id_relations_norm(isnan(id_relations_norm))=0;

fprintf('\n\n   & ');
for ii=1:7
    fprintf('%s & ',Materials_names{ii});
end
fprintf('\b\b\\\\ \\hline\n');
for ii=1:7
    fprintf(['    %s & ',repmat('%.3f & ',1,7),'\b\b \\\\ \\hline\n'],Materials_names{ii},id_relations_norm(ii,:))
end

%%
fprintf(['id_relations_norm = [',repmat(['[',repmat('%.4f, ',1,10),'\b\b];...\n'],1,10),'\b\b\b\b\b];\n'],id_relations_norm);
