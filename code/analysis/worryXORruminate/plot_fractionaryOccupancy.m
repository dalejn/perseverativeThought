%% Fractional occupancy

FractionalOccupancy_neutral = zeros(nsubjs,numClusters);
FractionalOccupancy_worry = zeros(nsubjs,numClusters);
FractionalOccupancy_ruminate = zeros(nsubjs,numClusters);
FractionalOccupancy_neutral_thought = zeros(nsubjs,numClusters);
FractionalOccupancy_neutral_oneBack = zeros(nsubjs,numClusters);
FractionalOccupancy_worry_thought = zeros(nsubjs,numClusters);
FractionalOccupancy_worry_oneBack = zeros(nsubjs,numClusters);
FractionalOccupancy_ruminate_thought = zeros(nsubjs,numClusters);
FractionalOccupancy_ruminate_oneBack = zeros(nsubjs,numClusters);

% scanInd = neutral or worry/ruminate
% nbackInd = thought cue or 1-back
for N = 1:nsubjs
    for K = 1:numClusters
        for J = 0:2
            subjPartition = partition(subjInd == N & scanInd == J);
            TwoBackBlock = nbackInd(subjInd == N & scanInd == J); % trial type for a subject and condition
            if J == 0
                FractionalOccupancy_neutral(N,K) = sum(subjPartition==K)/ nTR_all(N);
            elseif J==1
                FractionalOccupancy_worry(N,K) = sum(subjPartition==K)/ nTR_all(N);
            elseif J==2
                FractionalOccupancy_ruminate(N,K) = sum(subjPartition==K)/ nTR_all(N);
            end
            for I = 0:1
                if I == 0 && J == 0
                    FractionalOccupancy_neutral_thought(N,K) = sum(subjPartition(TwoBackBlock==I)==K) / length(subjPartition(TwoBackBlock==I));
                elseif I==1 && J==0
                    FractionalOccupancy_neutral_oneBack(N,K) = sum(subjPartition(TwoBackBlock==I)==K) / length(subjPartition(TwoBackBlock==I));
                elseif I==0 && J==1
                    FractionalOccupancy_worry_thought(N,K) = sum(subjPartition(TwoBackBlock==I)==K) / length(subjPartition(TwoBackBlock==I));
                elseif I==1 && J==1
                    FractionalOccupancy_worry_oneBack(N,K) = sum(subjPartition(TwoBackBlock==I)==K) / length(subjPartition(TwoBackBlock==I));
                elseif I==0 && J==2
                    FractionalOccupancy_ruminate_thought(N,K) = sum(subjPartition(TwoBackBlock==I)==K) / length(subjPartition(TwoBackBlock==I));
                elseif I==1 && J==2
                    FractionalOccupancy_ruminate_oneBack(N,K) = sum(subjPartition(TwoBackBlock==I)==K) / length(subjPartition(TwoBackBlock==I));
                end
            end
        end
    end
end
    
% save(fullfile(savedir,['FractionalOccupancy_k',num2str(numClusters),'.mat']),'FractionalOccupancy')
%% Plot fractional occupancy
savedir = '/Users/dalejn/Desktop/Dropbox/Projects/inProgress/2023-02-perseverantThought/results/';

for j = 1:3
    fig=figure(j);
    if j==1
        h = barh(1:nsubjs,FractionalOccupancy_neutral,'stacked');
        type = 'neutral';
        title(['Fractional occupancy during ' type])
    elseif j==2
        h = barh(1:nsubjs,FractionalOccupancy_worry,'stacked');
        type = 'worry';
        title(['Fractional occupancy during ' type])
    elseif j==3
        h = barh(1:nsubjs,FractionalOccupancy_ruminate,'stacked');
        type = 'ruminate';
        title(['Fractional occupancy during ' type])
    end
    for i=1:size(myColorMap,1)
        h(i).FaceColor = myColorMap(i,:);
    end
    colormap(myColorMap);
    set(gca,'XTick',[]); %which will get rid of all the markings for the y axis
    set(gca,'YTick',[]); %which will get rid of all the markings for the y axis
    set(gca,'Yticklabel',[]) 
    set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
    h = axes(fig,'visible','off'); 
    % c = lcolorbar(clusterNames);
    % set(c, 'Position', [0.92 0.168 0.022 0.7])
    saveas(fig,fullfile(basedir,'results/',['fractionalOccupancy','_k',num2str(numClusters),'_',type,'.pdf']));
end

%%
for condition = 1:2 % 1 for worry, 2 for ruminate
    types = {'neutral', 'worry', 'rumination'}
    %% plot fractional occupancy group differences (analysis i)
    
    
    % FractionalOccupancy_neutral = zeros(nsubjs,numClusters);
    % FractionalOccupancy_worry_ruminate = zeros(nsubjs,numClusters);
    % FractionalOccupancy_neutral_thought = zeros(nsubjs,numClusters);
    % FractionalOccupancy_neutral_oneBack = zeros(nsubjs,numClusters);
    % FractionalOccupancy_worry_ruminate_thought = zeros(nsubjs,numClusters);
    % FractionalOccupancy_worry_ruminate_oneBack = zeros(nsubjs,numClusters);
    
    thought_types = {'Neutral', char(types(condition+1))};
    group_types = {'control', 'clinical'};
    
    if condition==1
        FractionalOccupancy_worry_ruminate_thought = FractionalOccupancy_worry_thought;
        FractionalOccupancy_worry_ruminate_oneBack = FractionalOccupancy_worry_oneBack;
    elseif condition==2
        FractionalOccupancy_worry_ruminate_thought = FractionalOccupancy_ruminate_thought; 
        FractionalOccupancy_worry_ruminate_oneBack = FractionalOccupancy_ruminate_oneBack;
    end
    
    fig = figure(1)
    for J = 0
        for K=1:2
            thought_type = thought_types{K};
    
            if strcmpi(thought_type, 'Neutral')
                subplot(2,1,1)
                x = FractionalOccupancy_neutral_thought(groupInd==J,:);
                y = FractionalOccupancy_neutral_thought(groupInd==J+1,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters);
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters)
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                disp(1)
                [y,x] = find(p_adj<0.05)
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title("Fractional occupancy in neutral thought difference between clinical and control groups")

    %             for k = 1:numClusters
    %                 set(b(k), 'FaceColor', myColorMap(k,:));
    %             end
            else
                subplot(2,1,2)
                x = FractionalOccupancy_worry_ruminate_thought(groupInd==J,:);
                y = FractionalOccupancy_worry_ruminate_thought(groupInd==J+1,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters);
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters)
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                disp(2)
                [y,x] = find(p_adj<0.05)
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title("Fractional occupancy in worry/rumination difference between clinical and control groups")

    %             for k = 1:numClusters
    %                 set(b(k), 'FaceColor', myColorMap(k,:));
    %             end
            end
            
        end
    end
    colormap(myColorMap)
    h = axes(fig,'visible','off'); 
    h.Title.Visible = 'on';
    h.XLabel.Visible = 'on';
    h.YLabel.Visible = 'on';
    h.YLabel.Position(1) = -0.07;
    ylabel(h,'Clinical-control');
    xlabel(h,'Brain State','FontWeight','bold');
    
    saveas(fig,fullfile(savedir,strcat('fractionalOccupancy_', thought_type), '/',['between','diff_','fractionalOccupancy_k',num2str(numClusters),'.pdf']),'pdf');
    % see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
    % directed network
    
    %% plot fractional occupancy within group (analysis ii)
    thought_types = {'Neutral', char(types(condition+1))};
    group_types = {'control', 'clinical'};
    f = figure;
        
    for J = 0:1
        group_type = group_types{J+1};
        thought_type = thought_types{2};
    
        subplot(2,1,J+1)
        x = FractionalOccupancy_neutral_thought(groupInd==J,:);
        y = FractionalOccupancy_worry_ruminate_thought(groupInd==J,:);
        [~,p] = ttest2(y,x);
        find(p<0.05/numClusters);
        b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
        set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
        [y,x] = find(p<0.05/numClusters)
        text(x-.12,y-1,'*','Color','k');
        [~, ~, ~, p_adj] = fdr_bh(p)
        disp(3)
        [y,x] = find(p_adj<0.05)
        text(x-.12,y-1.01,'x','Color','k');
        b.CData = myColorMap;
        if J==0
            title(["Difference in fractional occupancy between neutral and " thought_type " within control group"])
        elseif J==1
            title(["Difference in fractional occupancy between neutral and " thought_type " within clinical group"])
        end
    end
    colormap(myColorMap)
    h = axes(f,'visible','off'); 
    h.Title.Visible = 'on';
    h.XLabel.Visible = 'on';
    h.YLabel.Visible = 'on';
    h.YLabel.Position(1) = -0.07;

    ylabel(h,[char(types(condition+1)) ' - neutral']);
    xlabel(h,'Brain State','FontWeight','bold');
    % c = lcolorbar(clusterNames);
    % set(c, 'Position', [0.92 0.168 0.022 0.7])
    saveas(f,fullfile(savedir, strcat('fractionalOccupancy_', thought_type), '/',['within','diff_', group_type,'_fractionalOccupancy_k',num2str(numClusters),'.pdf']),'pdf');
    
    % see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
    % directed network
    
    %% %% plot fractional occupancy group differences (analysis iii)
    
    thought_types = {'Neutral', char(types(condition+1))};
    group_types = {'control', 'clinical'};
    for J = 0
        for K=1:2
            thought_type = thought_types{K};
            fig = figure(K)
            if strcmpi(thought_type, 'Neutral')
                subplot(2,1,1)
                x = FractionalOccupancy_neutral_thought(groupInd==J,:);
                y = FractionalOccupancy_neutral_oneBack(groupInd==J,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters);
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters);
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                disp(4)
                [y,x] = find(p_adj<0.05)
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title(["Difference in fractional occupancy between 1-back and " thought_type " in control group"])                
    %             for k = 1:numClusters
    %                 set(b(k), 'FaceColor', myColorMap(k,:));
    %             end

                subplot(2,1,2)
                x = FractionalOccupancy_neutral_thought(groupInd==J+1,:);
                y = FractionalOccupancy_neutral_oneBack(groupInd==J+1,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters)
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters);
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                [y,x] = find(p_adj<0.05)
                disp(5)
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title(["Difference in fractional occupancy between 1-back and " thought_type " in clinical group"])                

    %             for k = 1:numClusters
    %                 set(b(k), 'FaceColor', myColorMap(k,:));
    %             end
                colormap(myColorMap)
                h = axes(fig,'visible','off'); 
                h.Title.Visible = 'on';
                h.XLabel.Visible = 'on';
                h.YLabel.Visible = 'on';
                h.YLabel.Position(1) = -0.07;
                ylabel(h,['1-back - ' char(types(condition+1))]);
                xlabel(h,'Brain State','FontWeight','bold');

                saveas(fig,fullfile(savedir, 'fractionalOccupancy_neutral/',['within',strcat('_diff_',thought_type,'ToOneBack'),'dwellTime_k',num2str(numClusters),'.pdf']),'pdf');
            else
                subplot(2,1,1)
                x = FractionalOccupancy_worry_ruminate_thought(groupInd==J,:);
                y = FractionalOccupancy_worry_ruminate_oneBack(groupInd==J,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters);
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters);
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                disp(6)
                [y,x] = find(p_adj<0.05)
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title(["Difference in fractional occupancy between 1-back and " thought_type " in control group"])                

                subplot(2,1,2)
                x = FractionalOccupancy_worry_ruminate_thought(groupInd==J+1,:);
                y = FractionalOccupancy_worry_ruminate_oneBack(groupInd==J+1,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters);
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters);
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                disp(7)
                [y,x] = find(p_adj<0.05)
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title(["Difference in fractional occupancy between 1-back and " thought_type " in clinical group"])

                colormap(myColorMap)
                h = axes(fig,'visible','off'); 
                h.Title.Visible = 'on';
                h.XLabel.Visible = 'on';
                h.YLabel.Visible = 'on';
                h.YLabel.Position(1) = -0.07;
                ylabel(h,['1-back - ' char(types(condition+1))]);
                xlabel(h,'Brain State','FontWeight','bold');
                % c = lcolorbar(clusterNames);
                % set(c, 'Position', [0.93 0.168 0.022 0.7])

                saveas(fig,fullfile(savedir, strcat('fractionalOccupancy_', thought_type), '/',['within', strcat('_diff_',thought_type,'ToOneBack'),'fractionalOccupancy_k',num2str(numClusters),'.pdf']),'pdf');
            end     
        end
    end
end