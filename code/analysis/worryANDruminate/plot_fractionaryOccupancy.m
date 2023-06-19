%% Fractional occupancy

FractionalOccupancy_neutral = zeros(nsubjs,numClusters);
FractionalOccupancy_worry_ruminate = zeros(nsubjs,numClusters);
FractionalOccupancy_neutral_thought = zeros(nsubjs,numClusters);
FractionalOccupancy_neutral_oneBack = zeros(nsubjs,numClusters);
FractionalOccupancy_worry_ruminate_thought = zeros(nsubjs,numClusters);
FractionalOccupancy_worry_ruminate_oneBack = zeros(nsubjs,numClusters);

% scanInd = neutral or worry/ruminate
% nbackInd = thought cue or 1-back
for N = 1:nsubjs
    for K = 1:numClusters
        for J = 0:1
            subjPartition = partition(subjInd == N & scanInd == J);
            TwoBackBlock = nbackInd(subjInd == N & scanInd == J);
            if J == 0
                FractionalOccupancy_neutral(N,K) = sum(subjPartition==K)/ nTR_all(N);
            else
                FractionalOccupancy_worry_ruminate(N,K) = sum(subjPartition==K)/ nTR_all(N);
            end
            for I = 0:1
                if I == 0 && J == 0
                    FractionalOccupancy_neutral_thought(N,K) = sum(subjPartition(TwoBackBlock==I)==K) / length(subjPartition(TwoBackBlock==J));
                elseif I==1 && J==0
                    FractionalOccupancy_neutral_oneBack(N,K) = sum(subjPartition(TwoBackBlock==I)==K) / length(subjPartition(TwoBackBlock==J));
                elseif I==0 && J==1
                    FractionalOccupancy_worry_ruminate_thought(N,K) = sum(subjPartition(TwoBackBlock==I)==K) / length(subjPartition(TwoBackBlock==J));
                elseif I==1 && J==1
                    FractionalOccupancy_worry_ruminate_oneBack(N,K) = sum(subjPartition(TwoBackBlock==I)==K) / length(subjPartition(TwoBackBlock==J));
                end
            end
        end
    end
end

% save(fullfile(savedir,['FractionalOccupancy_k',num2str(numClusters),'.mat']),'FractionalOccupancy')
%% Plot fractional occupancy
for j = 1:2
    fig=figure(j);
    if j==1
        h = barh(1:nsubjs,FractionalOccupancy_neutral,'stacked');
        type = 'neutral';
    elseif j==2
        h = barh(1:nsubjs,FractionalOccupancy_worry_ruminate,'stacked');
        type = 'worry_ruminate';
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
    c = lcolorbar(clusterNames);
    set(c, 'Position', [0.92 0.168 0.022 0.7])
    saveas(fig,fullfile(basedir,'results/',['fractionalOccupancy','_k',num2str(numClusters),'_',type,'.pdf']));
end

%% plot fractional occupancy group differences (analysis i)


% FractionalOccupancy_neutral = zeros(nsubjs,numClusters);
% FractionalOccupancy_worry_ruminate = zeros(nsubjs,numClusters);
% FractionalOccupancy_neutral_thought = zeros(nsubjs,numClusters);
% FractionalOccupancy_neutral_oneBack = zeros(nsubjs,numClusters);
% FractionalOccupancy_worry_ruminate_thought = zeros(nsubjs,numClusters);
% FractionalOccupancy_worry_ruminate_oneBack = zeros(nsubjs,numClusters);

thought_types = {'Neutral', 'Worry or rumination'};
group_types = {'control', 'clinical'};
fig = figure(1)
for J = 0
    for K=1:2
        thought_type = thought_types{K};

        if strcmpi(thought_type, 'Neutral')
            subplot(2,1,1)
            x = FractionalOccupancy_neutral_thought(groupInd==J,:);
            y = FractionalOccupancy_neutral_thought(groupInd==J+1,:);
            [~,p] = ttest2(y,x);
            find(p<0.05/numClusters)
            b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
            b.CData = myColorMap;
%             for k = 1:numClusters
%                 set(b(k), 'FaceColor', myColorMap(k,:));
%             end
        else
            subplot(2,1,2)
            x = FractionalOccupancy_worry_ruminate_thought(groupInd==J,:);
            y = FractionalOccupancy_worry_ruminate_thought(groupInd==J+1,:);
            [~,p] = ttest2(y,x);
            find(p<0.05/numClusters)
            b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
            b.CData = myColorMap;
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
ylabel(h,'Difference in Fractional Occupancy (clinical-control)','FontWeight','bold');
xlabel(h,'Brain State','FontWeight','bold');
c = lcolorbar(clusterNames);
set(c, 'Position', [0.93 0.168 0.022 0.7])

saveas(fig,fullfile(savedir,['between','diff_','fractionalOccupacy_k',num2str(numClusters),'.pdf']),'pdf');
% see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
% directed network

%% plot fractional occupancy within group (analysis ii)
thought_types = {'Neutral', 'Worry or rumination'};
group_types = {'control', 'clinical'};
f = figure;
    
for J = 0:1
    group_type = group_types{J+1};
    thought_type = thought_types{J+1};

    subplot(2,1,J+1)
    x = FractionalOccupancy_neutral_thought(groupInd==J,:);
    y = FractionalOccupancy_worry_ruminate_thought(groupInd==J,:);
    [~,p] = ttest2(y,x);
    find(p<0.05/numClusters)
    b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
    b.CData = myColorMap;

end
colormap(myColorMap)
h = axes(f,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';
ylabel(h,'Difference in Fractional Occupancy (worry or ruminate-neutral)','FontWeight','bold');
xlabel(h,'Brain State','FontWeight','bold');
c = lcolorbar(clusterNames);
set(c, 'Position', [0.92 0.168 0.022 0.7])
    saveas(f,fullfile(savedir,['within','diff_', group_type,'_fractionalOccupancy_k',num2str(numClusters),'.pdf']),'pdf');

% see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
% directed network

%% %% plot fractional occupancy group differences (analysis iii)

thought_types = {'Neutral', 'Worry or rumination'};
group_types = {'control', 'clinical'};
fig = figure(1)
for J = 0
    for K=1:2
        thought_type = thought_types{K};

        if strcmpi(thought_type, 'Neutral')
            subplot(2,1,1)
            x = FractionalOccupancy_worry_ruminate_thought(groupInd==J,:);
            y = FractionalOccupancy_worry_ruminate_oneBack(groupInd==J,:);
            [~,p] = ttest2(y,x);
            find(p<0.05/numClusters)
            b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
            b.CData = myColorMap;
%             for k = 1:numClusters
%                 set(b(k), 'FaceColor', myColorMap(k,:));
%             end
        else
            subplot(2,1,2)
            x = FractionalOccupancy_worry_ruminate_thought(groupInd==J+1,:);
            y = FractionalOccupancy_worry_ruminate_oneBack(groupInd==J+1,:);
            [~,p] = ttest2(y,x);
            find(p<0.05/numClusters)
            b=bar(nanmean(y,1)-nanmean(x,1), 'FaceColor', 'flat');
            b.CData = myColorMap;
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
ylabel(h,'Difference in Fractional Occupancy (1-back-worry/ruminate)','FontWeight','bold');
xlabel(h,'Brain State','FontWeight','bold');
c = lcolorbar(clusterNames);
set(c, 'Position', [0.93 0.168 0.022 0.7])

saveas(fig,fullfile(savedir,['within','worry_ruminate_diff_thoughtOneBack','fractionalOccupacy_k',num2str(numClusters),'.pdf']),'pdf');