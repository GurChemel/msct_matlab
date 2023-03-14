clear; clc;

files_path = 'C:\Users\gchem\Downloads\full_body\';
file_list = dir([files_path,'*.txt']);
full_body_raw = zeros(271,666,length(file_list));
% files_path = 'C:\Users\gchem\Downloads\female_body\';
% file_list = dir([files_path,'*.txt']);
% full_body_raw = zeros(293,650,length(file_list));

fprintf('File Parsing   0.00%%..\n');
for ii = 1:length(file_list)
    fprintf('\b\b\b\b\b\b\b\b\b\b%6.2f%%..\n',100*(ii/length(file_list)));
    full_body_raw(:,:,ii) = readmatrix([files_path,file_list(ii).name]);
end
fprintf('\b Finished.\n');
%%
id_to_comp = readmatrix([files_path,'..\id_to_comp.txt']);

hash_table = -1*ones(1,4302);
hash_table(id_to_comp(:,2)+1392) = id_to_comp(:,3);

full_body_id = hash_table(full_body_raw+1392);
full_body_id = uint8(full_body_id);

%%
figure;
subplot(2,1,1); imshow(squeeze(sum(full_body_raw,1)),[])
subplot(2,1,2); imshow(squeeze(sum(full_body_raw,2)),[])

%%
figure;
full_body_norm = flipud(squeeze(sum(full_body_raw,2))');
full_body_norm = full_body_norm - min(full_body_norm,[],'all');
full_body_norm = (full_body_norm/max(full_body_norm,[],'all'));
imshow(full_body_norm)
imwrite(full_body_norm,"ImagesForPresentations/full_body_sideways.png");
%%
xcat_materials = load('xcat_materials.mat');
Materials_names = {'Air','Adipose_tissue','Water','Muscle','dry_rib_water','Blood','Red_marrow'};
num_mats = length(Materials_names);

xcat_id = full_body_id+1;

relevant_ids = zeros(1,num_mats);
xcat_shorten_ids = -1*ones(size(xcat_id));
for ii = 1:num_mats
    relevant_ids(ii) = xcat_materials.(Materials_names{ii});
    xcat_shorten_ids(xcat_id==relevant_ids(ii))=ii;
end
%%
id_relations = zeros(num_mats);
fprintf('ID Parsing   0.00%%..\n');
cnt = 1;
printed_val = 0;
max_cnt = (size(xcat_shorten_ids,1)-2)*(size(xcat_shorten_ids,2)-2)*(size(xcat_shorten_ids,3)-2);
for row=2:(size(xcat_shorten_ids,1)-1)
    for col=2:(size(xcat_shorten_ids,2)-1)
        for slice=2:(size(xcat_shorten_ids,3)-1)
            print_val = floor(10000*(cnt/max_cnt))/100;
            if print_val~=printed_val
                fprintf('\b\b\b\b\b\b\b\b\b\b%6.2f%%..\n',print_val);
                printed_val = print_val;
            end
            cnt = cnt + 1;
            cur_id = xcat_shorten_ids(row,col,slice);
            if cur_id == -1
                continue
            end
            if xcat_shorten_ids(row-1,col,slice) ~= -1
                id_relations(cur_id,xcat_shorten_ids(row-1,col,slice)) = id_relations(cur_id,xcat_shorten_ids(row-1,col,slice)) + 1;
            end
            if xcat_shorten_ids(row+1,col,slice) ~= -1
                id_relations(cur_id,xcat_shorten_ids(row+1,col,slice)) = id_relations(cur_id,xcat_shorten_ids(row+1,col,slice)) + 1;
            end
            if xcat_shorten_ids(row,col-1,slice) ~= -1
                id_relations(cur_id,xcat_shorten_ids(row,col-1,slice)) = id_relations(cur_id,xcat_shorten_ids(row,col-1,slice)) + 1;
            end
            if xcat_shorten_ids(row,col+1,slice) ~= -1
                id_relations(cur_id,xcat_shorten_ids(row,col+1,slice)) = id_relations(cur_id,xcat_shorten_ids(row,col+1,slice)) + 1;
            end
            if xcat_shorten_ids(row,col,slice-1) ~= -1
                id_relations(cur_id,xcat_shorten_ids(row,col,slice-1)) = id_relations(cur_id,xcat_shorten_ids(row,col,slice-1)) + 1;
            end
            if xcat_shorten_ids(row,col,slice+1) ~= -1
                id_relations(cur_id,xcat_shorten_ids(row,col,slice+1)) = id_relations(cur_id,xcat_shorten_ids(row,col,slice+1)) + 1;
            end
        end
    end
end
fprintf('\b Finished.\n');
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
clear; clc;

load('full_male_body_id_relations.mat');
full_male_body_id_relations = id_relations_norm;
load('full_female_body_id_relations.mat');
full_female_body_id_relations = id_relations_norm;
clear id_relations_norm;

percentage_helper = (full_male_body_id_relations+full_female_body_id_relations)/2;
percentage_helper(percentage_helper<0.001)=Inf;

figure;
Materials_names = {'Air','Adipose Tissue','Water','Muscle','Bone','Blood','Red Marrow'};
plot(100*(full_male_body_id_relations-full_female_body_id_relations)./percentage_helper)
legend(Materials_names)
% for ii=1:7
%     subplot(2,4,ii);
%     plot([full_male_body_id_relations(:,ii)-full_female_body_id_relations(:,ii)]);
%     title(Materials_names{ii});
% end