clear; clc;
date_str = '220325';
load(['OptimResults_',date_str,'.mat']);
titles_vec = {'Loss','Epsilon','Delta'};
legend_vec = {'Init From Direct Only       - Multispectral      Optimization',...
              'Init From Direct And Global - Multispectral      Optimization',...
              'Init From Direct Only       - Energy Accumulated Optimization',...
              'Init From Direct And Global - Energy Accumulated Optimization'};

%%
plot_0_or_save_1 = 1;

if plot_0_or_save_1==0
    figure;
    for ii=0:2
        subplot(2,2,ii+1);
        plot(OptimResults(:,(ii*4)+(1:4)),'LineWidth',3);
        title(titles_vec{ii+1})
    end
    % Legend:
    colors = [[0 0.4470 0.7410];...	'#0072BD'	
              [0.8500 0.3250 0.0980];...	'#D95319'	
              [0.9290 0.6940 0.1250];...	'#EDB120'	
              [0.4940 0.1840 0.5560]]; %	'#7E2F8E'	
          
    reorder_vec = [1,3,2,4];
    subplot(2,2,4)
    for ii=1:4
        plot(zeros(1,4),zeros(1,4),'LineWidth',3,'Color',colors(reorder_vec(ii),:)); hold on;
        legend_vec_reordered{ii} = legend_vec{reorder_vec(ii)};
    end
    legend(legend_vec_reordered,'FontName','FixedWidth','FontSize', 10)
end

if plot_0_or_save_1==1
    for ii=1:2
        h1 = figure;
        plot(OptimResults(:,(ii*4)+[2,4]),'LineWidth',3);
        set(gca,'FontSize',15);
        xt = xticks;
        xticks(xt(1:2:end));
        yt = yticks;
        yticks(yt(1:2:end));
        xlabel('Iteration')
        if ii==1
            ylabel('\epsilon_{atten}');
            legend('Multi Spectral','Energy Accumulated');
        else
            ylabel('\delta_{atten}');
            legend('Multi Spectral','Energy Accumulated','Location','southeast');
        end
        set(gca,'FontSize',15);
%         xt = xticks;
%         xticks(xt(1:2:end));
        yt = yticks;
        yticks(yt(1:2:end));
        grid on;
%         saveas(h1,sprintf('ImagesForPresentations/220824_%s_Optim_type_%s.png',date_str,titles_vec{ii+1}))
%         close(h1)
    end
end

%%

% if plot_0_or_save_1==1
    for ii=1:2
        h1 = figure;
        plot(OptimResults(:,(ii*4)+[2,4]),'LineWidth',3);
        set(gca,'FontSize',15);
        xt = xticks;
        xticks(xt(1:2:end));
        yt = yticks;
        yticks(yt(1:2:end));
        xlabel('Iteration')
        if ii==1
            ylabel('\epsilon_{atten}');
%             legend('Multi Spectral','Energy Accumulated');
        else
            ylabel('\delta_{atten}');
%             legend('Multi Spectral','Energy Accumulated','Location','southeast');
        end
%         set(gca,'FontSize',20);
%         xt = xticks;
%         xticks(xt(1:2:end));
        yt = yticks;
        yticks(yt(1:2:end));
        grid on;
%         saveas(h1,sprintf('ImagesForPresentations/220824_%s_Optim_type_%s.png',date_str,titles_vec{ii+1}))
%         close(h1)
    end
% end
%%
plot_0_or_save_1 = 0;

color_map = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980]];
color_map = [color_map ; brighten(color_map,.8)];
    
if plot_0_or_save_1==0
    figure;
    for ii=0:2
        subplot(2,2,ii+1);
        for jj=1:4
            plot(OptimResults(:,(ii*4)+jj),'LineWidth',3,'Color',color_map(jj,:)); hold on;
        end
        title(titles_vec{ii+1})
    end
    
    % Legend:
    reorder_vec = [1,3,2,4];
    subplot(2,2,4)
    for ii=1:4
        plot(zeros(1,4),zeros(1,4),'LineWidth',3,'Color',color_map(reorder_vec(ii),:)); hold on;
        legend_vec_reordered{ii} = legend_vec{reorder_vec(ii)};
    end
    legend(legend_vec_reordered,'FontName','FixedWidth','FontSize', 10)
end

if plot_0_or_save_1==1
%     for ii=0:2
%         h1 = figure;
%         for jj=1:4
%             plot(OptimResults(:,(ii*4)+jj),'LineWidth',3,'Color',color_map(jj,:)); hold on;
%         end
%         xlim([1,size(OptimResults,1)])
%         set(gca,'FontSize',15);
%         xt = xticks;
%         xticks(xt(1:2:end));
%         yt = yticks;
%         yticks(yt(1:2:end));
%         grid on;
%         saveas(h1,sprintf('ImagesForPresentations/%s_Optim_type_%s.png',date_str,titles_vec{ii+1}))
%         close(h1)
%     end
    
    h1 = figure;
    for ii=1:2
        subplot(1,2,ii);
        for jj=[3,4,1,2] %1:4
            plot(OptimResults(:,(ii*4)+jj),'LineWidth',3,'Color',color_map(jj,:)); hold on;
        end
        xlim([1,size(OptimResults,1)])
        xlabel('Iteration')
        if ii==1
            ylabel('\epsilon');
            legend('Scenario 1','Scenario 2','Scenario 3','Scenario 4');
        else
            ylabel('\delta_{mass}');
            legend('Scenario 1','Scenario 2','Scenario 3','Scenario 4','Location','southeast');
        end
        set(gca,'FontSize',15);
%         xt = xticks;
%         xticks(xt(1:2:end));
        yt = yticks;
        yticks(yt(1:2:end));
        grid on;
    end
    saveas(h1,sprintf('ImagesForPresentations/%s_Optim_type_%s_and_%s.png',date_str,titles_vec{2},titles_vec{3}))
%     close(h1)

%     reorder_vec = [1,3,2,4];
%     h1 = figure;
%     for ii=1:4
%         plot(zeros(1,4),zeros(1,4),'LineWidth',3,'Color',color_map(reorder_vec(ii),:)); hold on;
%         legend_vec_reordered{ii} = legend_vec{reorder_vec(ii)};
%     end
%     legend(legend_vec_reordered,'FontName','FixedWidth','FontSize', 14)
%     saveas(h1,sprintf('ImagesForPresentations/%s_Optim_type_Legend.png',date_str))
%     close(h1)
end