%% Dwell time

dwellTime_neutral = zeros(nsubjs,numClusters);
dwellTime_worry = zeros(nsubjs,numClusters);
dwellTime_ruminate = zeros(nsubjs,numClusters);
dwellTime_neutral_thought = zeros(nsubjs,numClusters);
dwellTime_neutral_oneBack = zeros(nsubjs,numClusters);
dwellTime_worry_thought = zeros(nsubjs,numClusters);
dwellTime_worry_oneBack = zeros(nsubjs,numClusters);
dwellTime_ruminate_thought = zeros(nsubjs,numClusters);
dwellTime_ruminate_oneBack = zeros(nsubjs,numClusters);
% calculate dwell time, i.e. average number of subsequent TRs that each state lasts for
% store both mean and median, because dwell time may not be normally distributed
TR = 3;     % TR length
DwellTimeMean = zeros(nsubjs,numClusters);
% DwellTimeMedian = zeros(nsubjs,numClusters);
% RunRate = zeros(nsubjs,numClusters);

% scanInd = neutral or worry/ruminate
% nbackInd = thought cue or 1-back
for N = 1:nsubjs
    for J = 0:2
        subjPartition = partition(subjInd == N & scanInd == J);
        TwoBackBlock = nbackInd(subjInd == N & scanInd == J);
        if J == 0
            dwellTime_neutral(N,:) = CALC_DWELL_TIME(subjPartition,numClusters);
        elseif J == 1
            dwellTime_worry(N,:) = CALC_DWELL_TIME(subjPartition,numClusters);
        elseif J == 2
            dwellTime_ruminate(N,:) = CALC_DWELL_TIME(subjPartition,numClusters);
        end
        for I = 0:1
            if I == 0 && J == 0
                dwellTime_neutral_thought(N,:) = CALC_DWELL_TIME(subjPartition(TwoBackBlock==I),numClusters);
            elseif I==1 && J==0
                dwellTime_neutral_oneBack(N,:) = CALC_DWELL_TIME(subjPartition(TwoBackBlock==I),numClusters);
            elseif I==0 && J==1
                dwellTime_worry_thought(N,:) = CALC_DWELL_TIME(subjPartition(TwoBackBlock==I),numClusters);
            elseif I==1 && J==1
                dwellTime_worry_oneBack(N,:) = CALC_DWELL_TIME(subjPartition(TwoBackBlock==I),numClusters);
            elseif I==0 && J==2
                dwellTime_ruminate_thought(N,:) = CALC_DWELL_TIME(subjPartition(TwoBackBlock==I),numClusters);
            elseif I==1 && J==2
                dwellTime_ruminate_oneBack(N,:) = CALC_DWELL_TIME(subjPartition(TwoBackBlock==I),numClusters);
            end
        end
    end
end

% for N = 1:nsubjs
%     [dt_mean,dt_median,~,~,n_runs] = CALC_DWELL_TIME(partition(subjInd' == N & scanInd == (i-1)),numClusters);
%     DwellTimeMean(N,:) = dt_mean*TR;        % store dwell time in seconds
%     DwellTimeMedian(N,:) = dt_median*TR;        % store dwell time in seconds
%     % store rate of appearance of runs, i.e. how many DMN runs appear in 1 minute
%     % first calculate runs/TR,then runs/sec, then runs/min
%     RunRate(N,:) = 60*(1/TR)*(n_runs/numTRs(i));   
% end

% save(fullfile(savedir,[scanlab{i},'DwellTime_k',num2str(numClusters),name_root,'.mat']),'DwellTimeMean','DwellTimeMedian','RunRate')
savedir = '/Users/dalejn/Desktop/Dropbox/Projects/inProgress/2023-02-perseverantThought/results/';

hexMap = GET_CLUSTER_COLORS(6);
for k = 1 : length(hexMap)
	thisCell = hexMap{k}
	r = hex2dec(thisCell(1:2))
	g = hex2dec(thisCell(3:4))
	b = hex2dec(thisCell(5:6))
	myColorMap(k, :) = [r, g, b]
end
myColorMap = myColorMap / 255; % Normalize to range 0-1


%% Plot dwell time
for j = 1:3
    fig=figure(j);
    if j==1
        h = barh(1:nsubjs,dwellTime_neutral,'stacked');
        type = 'neutral';
        title(['Dwell times during ' type])
    elseif j==2
        h = barh(1:nsubjs,dwellTime_worry,'stacked');
        type = 'worry';
        title(['Dwell times during ' type])
    elseif j==3
        h = barh(1:nsubjs,dwellTime_ruminate,'stacked');
        type = 'ruminate';
        title(['Dwell times during ' type])
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
    saveas(fig,fullfile(basedir,'results/',['dwellTime','_k',num2str(numClusters),'_',type,'.pdf']));
end

%% for worry and ruminate groups separately

for condition = 1:2

    %% plot dwell time group differences (analysis i)
    
    thought_types = {'Neutral', char(types(condition+1))};
    group_types = {'control', 'clinical'};
    
    if condition==1
        dwellTime_worry_ruminate_thought = dwellTime_worry_thought;
        dwellTime_worry_ruminate_oneBack = dwellTime_worry_oneBack;
    elseif condition==2
        dwellTime_worry_ruminate_thought = dwellTime_ruminate_thought; 
        dwellTime_worry_ruminate_oneBack = dwellTime_ruminate_oneBack;
    end
    
    fig = figure(1)
    for J = 0
        for K=1:2
            thought_type = thought_types{K};
    
            if strcmpi(thought_type, 'Neutral')
                subplot(2,1,1)
                x = dwellTime_neutral_thought(groupInd==J,:);
                y = dwellTime_neutral_thought(groupInd==J+1,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters)
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters)
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                [y,x] = find(p_adj<0.05);
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title("Dwell time in neutral thought difference between clinical and control groups")
    %             for k = 1:numClusters
    %                 set(b(k), 'FaceColor', myColorMap(k,:));
    %             end
            else
                subplot(2,1,2)
                x = dwellTime_worry_ruminate_thought(groupInd==J,:);
                y = dwellTime_worry_ruminate_thought(groupInd==J+1,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters)
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters)
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                [y,x] = find(p_adj<0.05);
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title(strcat("Dwell time in ",thought_type," difference between clinical and control groups"))
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
    ylabel(h,'Clinical-control','FontWeight','bold');
    xlabel(h,'Brain State','FontWeight','bold');
    
    saveas(fig,fullfile(savedir, strcat('dwellTime_', thought_type), '/',['between','diff_','dwellTime_k',num2str(numClusters),'.pdf']),'pdf');
    % see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
    % directed network
    
    %% plot dwell time within group (analysis ii)
    
    thought_types = {'Neutral', char(types(condition+1))};
    group_types = {'control', 'clinical'};
    f = figure;
        
    for J = 0:1
        group_type = group_types{J+1};
        thought_type = thought_types{2};
    
        subplot(2,1,J+1)
        x = dwellTime_neutral_thought(groupInd==J,:);
        y = dwellTime_worry_ruminate_thought(groupInd==J,:);
        [~,p] = ttest2(y,x);
        find(p<0.05/numClusters)
        b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
        set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
        [y,x] = find(p<0.05/numClusters)
        text(x-.12,y-1,'*','Color','k');
        [~, ~, ~, p_adj] = fdr_bh(p)
        [y,x] = find(p_adj<0.05);
        text(x-.12,y-1.01,'x','Color','k');
        b.CData = myColorMap;
        if J==0
            title(["Difference in dwell time between neutral and " thought_type " within control group"])
        elseif J==1
            title(["Difference in dwell time between neutral and " thought_type " within clinical group"])
        end
    end
    colormap(myColorMap)
    h = axes(f,'visible','off'); 
    h.Title.Visible = 'on';
    h.XLabel.Visible = 'on';
    h.YLabel.Visible = 'on';
    h.YLabel.Position(1) = -0.07;
    ylabel(h,[thought_type ' - neutral'],'FontWeight','bold');
    xlabel(h,'Brain State','FontWeight','bold');
    
    saveas(f,fullfile(savedir, strcat('dwellTime_', thought_type), '/',['within','diff_', group_type,'_dwellTime_k',num2str(numClusters),'.pdf']),'pdf');
    
    % see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
    % directed network
    
    %% %% plot dwell time group differences (analysis iii)
    
    thought_types = {'Neutral', char(types(condition+1))};
    group_types = {'control', 'clinical'};
    for J = 0
        for K=1:2
            thought_type = thought_types{K};
            fig = figure(K)
            if strcmpi(thought_type, 'Neutral')
                subplot(2,1,1)
                x = dwellTime_neutral_thought(groupInd==J,:);
                y = dwellTime_neutral_oneBack(groupInd==J,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters)
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters)
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                [y,x] = find(p_adj<0.05);
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title(["Difference in dwell time between 1-back and " thought_type " in control group"])
    %             for k = 1:numClusters
    %                 set(b(k), 'FaceColor', myColorMap(k,:));
    %             end

                subplot(2,1,2)
                x = dwellTime_neutral_thought(groupInd==J+1,:);
                y = dwellTime_neutral_oneBack(groupInd==J+1,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters)
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters)
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                [y,x] = find(p_adj<0.05);
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title(["Difference in dwell time between 1-back and " thought_type " in clinical group"])

                colormap(myColorMap)
                h = axes(fig,'visible','off'); 
                h.Title.Visible = 'on';
                h.XLabel.Visible = 'on';
                h.YLabel.Visible = 'on';
                h.YLabel.Position(1) = -0.07;
                ylabel(h,['1-back - ' thought_type],'FontWeight','bold');
                xlabel(h,'Brain State','FontWeight','bold');
                
                saveas(fig,fullfile(savedir, 'dwellTime_neutral/',['within',strcat('_diff_',thought_type,'ToOneBack'),'dwellTime_k',num2str(numClusters),'.pdf']),'pdf');
            else
                subplot(2,1,1)
                x = dwellTime_worry_ruminate_thought(groupInd==J,:);
                y = dwellTime_worry_ruminate_oneBack(groupInd==J,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters)
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters)
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                [y,x] = find(p_adj<0.05);
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title(["Difference in dwell time between 1-back and " thought_type " in control group"])                

                subplot(2,1,2)
                x = dwellTime_worry_ruminate_thought(groupInd==J+1,:);
                y = dwellTime_worry_ruminate_oneBack(groupInd==J+1,:);
                [~,p] = ttest2(y,x);
                find(p<0.05/numClusters)
                b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
                set(gca, 'XTick', 1:numClusters,'XTickLabel',clusterNames);
                [y,x] = find(p<0.05/numClusters)
                text(x-.12,y-1,'*','Color','k');
                [~, ~, ~, p_adj] = fdr_bh(p)
                [y,x] = find(p_adj<0.05);
                text(x-.12,y-1.01,'x','Color','k');
                b.CData = myColorMap;
                title(["Difference in dwell time between 1-back and " thought_type " in clinical group"])

                colormap(myColorMap)
                h = axes(fig,'visible','off'); 
                h.Title.Visible = 'on';
                h.XLabel.Visible = 'on';
                h.YLabel.Visible = 'on';
                h.YLabel.Position(1) = -0.07;
                ylabel(h,['1-back - ' thought_type],'FontWeight','bold');
                xlabel(h,'Brain State','FontWeight','bold');
                
                saveas(fig,fullfile(savedir, strcat('dwellTime_', thought_type), '/',['within', strcat('_diff_',thought_type,'ToOneBack'),'dwellTime_k',num2str(numClusters),'.pdf']),'pdf');
            end
        end
    end
end
