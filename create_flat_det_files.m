clear; clc;

src_loc = csvread('source_and_det_2\src_loc.csv');
cut_frames = csvread('source_and_det_2\cut_frames_2.csv');


det_loc = -src_loc;
v = repmat([0,0,1],180,1);
u = cross(det_loc,v);
u = u./repmat(vecnorm(u')',1,3);

%% Check det vs src

for ii=1:180
    if (abs(det_loc(ii,:)*u(ii,:)') > 1e-8)
        error('Dot product not zero for idx u(%d)',ii);
    end
    if (abs(det_loc(ii,:)*v(ii,:)') > 1e-8)
        error('Dot product not zero for idx v(%d)',ii);
    end
    if (abs(u(ii,:)*v(ii,:)') > 1e-8)
        error('Dot product not zero for idx u*v(%d)',ii);
    end
end

%%
figure;
for ii=1:180
    scatter3(src_loc(ii,1),src_loc(ii,2),src_loc(ii,3),'.'); hold on;
    scatter3(det_loc(ii,1),det_loc(ii,2),det_loc(ii,3),'.'); hold on;
    xlim([-280,280]);
    ylim([-280,280]);
    zlim([-2,2]);
    
    view(60,30);
    pause(0.1);
end

%%
% src_id = cut_frames(1:180,1);
% det_id_borders = cut_frames(1:180,[2,3])+1;
% det_angle_per_source = zeros(length(src_id),401);
% for loop_num=1:length(src_id)
%     if diff(det_id_borders(loop_num,:))>0
%         det_angle_per_source(loop_num,:) = det_id_borders(loop_num,1):det_id_borders(loop_num,2);
%     else
%         det_angle_per_source(loop_num,:) = [det_id_borders(loop_num,1):1800,1:det_id_borders(loop_num,2)];
%     end
% end
% 
% det_angle_per_source = det_angle_per_source-1;
% det_id_per_source = kron(det_angle_per_source,ones(1,80))*80+repmat(0:79',180,401);
% 
% 
