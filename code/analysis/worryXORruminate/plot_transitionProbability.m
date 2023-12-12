%% transition probabilities

% transition probability matrices:
% restTransitionProbabilityMats: each element reflects probability that
% state j occurs at the TR after state i, given that you are in state i (and not in the last TR of a scan)
% restTransitionProbabilityMatsNoPersist: each element reflectts
% probability that state j follows state i in the next state change,
% given that you are in state i (and not in the last TR of a scan). This
% *NoPersist version is used in the manuscript but both are interesting.

% get neutral thought transition probability matrices -- 2D is flattened by row for
% regressions

% scanInd = neutral or worry/ruminate
% nbackInd = thought cue or 1-back

[restTransitionProbability2D,restTransitionProbabilityMats] = GET_TRANS_PROBS(partition(scanInd == 0),subjInd(scanInd == 0));

neutralTransitionProbabilityNoPersist2D_oneBack = zeros(nsubjs,numClusters^2);
neutralTransitionProbabilityMatsNoPersist_oneBack = zeros(nsubjs,numClusters,numClusters);	% preallocate transition probability matrices
neutralTransitionProbabilityNoPersist2D_thought_cue = zeros(nsubjs,numClusters^2);
neutralTransitionProbabilityMatsNoPersist_thought_cue = zeros(nsubjs,numClusters,numClusters);	% preallocate transition probability matrices


for N = 1:nsubjs
	subjPartition = partition(subjInd == N & scanInd == 0);
    TwoBackBlock = nbackInd(subjInd == N & scanInd == 0);
    NotTwoBackBlock = ~nbackInd(subjInd == N & scanInd ==0);
    % add for loop for scanInd
	[neutralTransitionProbabilityMatsNoPersist_oneBack(N,:,:),neutralTransitionProbabilityNoPersist2D_oneBack(N,:)] = GET_BLOCK_TRANS_PROBS_NO_PERSIST(subjPartition,TwoBackBlock,numClusters);
    [neutralTransitionProbabilityMatsNoPersist_thought_cue(N,:,:),neutralTransitionProbabilityNoPersist2D_thought_cue(N,:)] = GET_BLOCK_TRANS_PROBS_NO_PERSIST(subjPartition,NotTwoBackBlock,numClusters);

end

%% get transition probabilities within two back blocks
% load data/HCP_TwoBackBlock.mat % load indicator vector specifying location of two back blocks

for condition = 1:2 % 1 for worry, 2 for ruminate
    types = {'Neutral', 'Worry', 'Rumination'}
    
    RWTransitionProbabilityNoPersist2D_oneBack = zeros(nsubjs,numClusters^2);
    RWTransitionProbabilityMatsNoPersist_oneBack = zeros(nsubjs,numClusters,numClusters);	% preallocate transition probability matrices
    RWTransitionProbabilityNoPersist2D_thought_cue = zeros(nsubjs,numClusters,numClusters);
    RWTransitionProbabilityMatsNoPersist_thought_cue = zeros(nsubjs,numClusters,numClusters);
    for N = 1:nsubjs
	    subjPartition = partition(subjInd == N & scanInd == condition);
        TwoBackBlock = nbackInd(subjInd == N & scanInd == condition);
        NotTwoBackBlock = ~nbackInd(subjInd == N & scanInd ==condition);
	    [RWTransitionProbabilityMatsNoPersist_oneBack(N,:,:),RWTransitionProbabilityNoPersist2D_oneBack(N,:)] = GET_BLOCK_TRANS_PROBS_NO_PERSIST(subjPartition,TwoBackBlock,numClusters);
        [RWTransitionProbabilityMatsNoPersist_thought_cue(N,:,:),RWTransitionProbabilityNoPersist2D_thought_cue(N,:)] = GET_BLOCK_TRANS_PROBS_NO_PERSIST(subjPartition,NotTwoBackBlock,numClusters);
    end
    grpAvgRest = squeeze(nanmean(neutralTransitionProbabilityMatsNoPersist_thought_cue,1)) .* ~eye(numClusters);
    % nans occur for transitions from states that are not present at all for a subject or within the tested blocks. 
    % transitions to that state are 0
    grpAvg2Back = squeeze(nanmean(RWTransitionProbabilityMatsNoPersist_oneBack,1)) .* ~eye(numClusters); 
    
    %% permutation testing to compare transition probability matrices
    disp('start permutation testing')
    nperms = 100000;
    pvals_twotail = PERM_TEST(RWTransitionProbabilityMatsNoPersist_thought_cue,neutralTransitionProbabilityMatsNoPersist_oneBack,nperms);
    
    %% plot transition probabilities group average
    maxVal = max(max([grpAvgRest,grpAvg2Back])); % sync color scales
    
    f = figure;
    
    subplot(1,3,1);
    imagesc(grpAvgRest);
    xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
    xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
    COLOR_TICK_LABELS(true,true,numClusters);
    ylabel('Current State'); xlabel('Next New State');
    title('Neutral');
    set(gca,'FontSize',8);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','arial');
    caxis([0 maxVal]); colorbar
    
    subplot(1,3,2);
    imagesc(grpAvg2Back);
    xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
    xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
    COLOR_TICK_LABELS(true,true,numClusters);
    
    ylabel('Current State'); xlabel('Next New State');
    title(char(types(condition+1)));
    set(gca,'FontSize',8);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','arial');
    caxis([0 maxVal]); colorbar
    
    subplot(1,3,3);
    nBackMinusRestTP = (grpAvg2Back-grpAvgRest);
    imagesc(nBackMinusRestTP.*~eye(numClusters)); colormap('plasma');
    xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
    yticks(1:numClusters); yticklabels(clusterNames); axis square
    ylabel('Current State'); xlabel('Next New State');
    
    sig_thresh = 0.05 / ((numClusters^2)-numClusters);      % bonferroni correction, for two-tailed p-values so only
    [y,x] = find(pvals_twotail.*~eye(numClusters) < sig_thresh);
    text(x-.12,y+.12,'*','Color','w');
    
    [~, ~, ~, adj_p] = fdr_bh(pvals_twotail.*~eye(numClusters));
    [y,x] = find(adj_p < 0.05);
    text(x-.12,y+.12,'x','Color','w');
    
    caxis_bound = max(max(abs(nBackMinusRestTP.*~eye(numClusters))));
    h = colorbar; ylabel(h, [char(types(condition+1)) ' - Neutral']); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
    COLOR_TICK_LABELS(true,true,numClusters);
    title([char(types(condition+1)) ' > Neutral']);
    set(gca,'FontSize',8);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','arial');
    
    f.PaperUnits = 'inches';
    f.PaperSize = [8 4];
    f.PaperPosition = [0 0 8 4];
    % plot transition probability matrices
    %saveas(f,fullfile(savedir,['Rest2BackTransProbs_k',num2str(numClusters),'.pdf']),'pdf');
    
    % see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
    % directed network
    
    %% plot transition probabilities group differences (analysis i)
    
    
    thought_types = {'Neutral', char(types(condition+1))};
    group_types = {'control', 'clinical'};
    for J = 0
        group_type = group_types{J+1};
        for K = 1:2
            thought_type = thought_types{K};
            
            if strcmpi(thought_type, 'Neutral')
                x = neutralTransitionProbabilityMatsNoPersist_thought_cue(groupInd==J,:);
                y = neutralTransitionProbabilityMatsNoPersist_thought_cue(groupInd==J+1,:);
                [~,p] = ttest2(y,x);
                x = reshape(nanmean(x,1), 6, 6);
                y = reshape(nanmean(y,1), 6, 6);
                maxVal = max(max([x,y])); % sync color scales
            else
                x = RWTransitionProbabilityMatsNoPersist_thought_cue(groupInd==J,:);
                y = RWTransitionProbabilityMatsNoPersist_thought_cue(groupInd==J+1,:);
                [~,p] = ttest2(y,x);
                x = reshape(nanmean(x,1), 6, 6);
                y = reshape(nanmean(y,1), 6, 6);
                maxVal = max(max([x,y])); % sync color scales
            end
            f = figure;
    
            subplot(1,3,1);
            imagesc(x);
            xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
            xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
            COLOR_TICK_LABELS(true,true,numClusters);
            ylabel('Current State'); xlabel('Next New State');
            title([thought_type, ', ', group_types{1}]);
            set(gca,'FontSize',8);
            set(gca,'TickLength',[0 0]);
            set(gca,'Fontname','arial');
            caxis([0 maxVal]); colorbar
    
            subplot(1,3,2);
            imagesc(y);
            xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
            xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
            COLOR_TICK_LABELS(true,true,numClusters);
    
            ylabel('Current State'); xlabel('Next New State');
            title([thought_type, ', ', group_types{2}]);
            set(gca,'FontSize',8);
            set(gca,'TickLength',[0 0]);
            set(gca,'Fontname','arial');
            caxis([0 maxVal]); colorbar
    
            subplot(1,3,3);
            nBackMinusRestTP = (y-x);
            imagesc(nBackMinusRestTP.*~eye(numClusters)); colormap('plasma');
            xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
            yticks(1:numClusters); yticklabels(clusterNames); axis square
            ylabel('Current State'); xlabel('Next New State');
            
            sig_thresh = 0.05 / ((numClusters^2)-numClusters);      % bonferroni correction, for two-tailed p-values so only
            [y,x] = find(reshape(p, 6, 6).*~eye(numClusters) < sig_thresh);
            text(x-.12,y+.12,'*','Color','w');
            
            [~, ~, ~, adj_p] = fdr_bh(reshape(p, 6, 6).*~eye(numClusters));
            [y,x] = find(adj_p < 0.05);
            text(x-.12,y+.12,'x','Color','w');
            
            caxis_bound = max(max(abs(nBackMinusRestTP.*~eye(numClusters))));
            h = colorbar; ylabel(h,[thought_type group_types{2} '-' thought_type group_types{1}]); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
            COLOR_TICK_LABELS(true,true,numClusters);
            title([thought_type, ', ', group_types{2}, '-', thought_type, ', ', group_types{1}]);
            set(gca,'FontSize',8);
            set(gca,'TickLength',[0 0]);
            set(gca,'Fontname','arial');

            if K==1
                sgtitle("Transition probabilities for neutral thought between clinical and control groups")
            elseif K==2 && condition==1
                sgtitle("Transition probabilities for worry between clinical and control groups")
            elseif K==2 && condition==2
                sgtitle("Transition probabilities for rumination between clinical and control groups")
            end

            f.PaperUnits = 'inches';
            f.PaperSize = [8 4];
            f.PaperPosition = [0 0 8 4];
            % plot transition probability matrices
            saveas(f,fullfile(savedir,strcat('transitionProbability_', thought_type), '/',['between_','diff_', regexprep(lower(thought_type), ' ', '_'),'_TransProbs_k',num2str(numClusters),'.pdf']),'pdf');
        end
    end
    % see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
    % directed network
    
    %% plot transition probabilities within group (analysis ii)
    thought_types = {'Neutral', char(types(condition+1))};
    group_types = {'control', 'clinical'};
    for J = 0:1
        group_type = group_types{J+1};
        thought_type = thought_types{J+1};
        f = figure;
        
        x = neutralTransitionProbabilityMatsNoPersist_thought_cue(groupInd==J,:);
        y = RWTransitionProbabilityMatsNoPersist_thought_cue(groupInd==J,:);
        [~,p] = ttest2(y,x);
        x = reshape(nanmean(x,1), 6, 6);
        y = reshape(nanmean(y,1), 6, 6);
        maxVal = max(max([x,y])); % sync color scales
    
        subplot(1,3,1);
        imagesc(x);
        xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
        xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
        COLOR_TICK_LABELS(true,true,numClusters);
        ylabel('Current State'); xlabel('Next New State');
        title(thought_types{1});
        set(gca,'FontSize',8);
        set(gca,'TickLength',[0 0]);
        set(gca,'Fontname','arial');
        caxis([0 maxVal]); colorbar
    
        subplot(1,3,2);
        imagesc(y);
        xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
        xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
        COLOR_TICK_LABELS(true,true,numClusters);
    
        ylabel('Current State'); xlabel('Next New State');
        title(thought_types{2});
        set(gca,'FontSize',8);
        set(gca,'TickLength',[0 0]);
        set(gca,'Fontname','arial');
        caxis([0 maxVal]); colorbar
    
        subplot(1,3,3);
        nBackMinusRestTP = (y-x);
        imagesc(nBackMinusRestTP.*~eye(numClusters)); colormap('plasma');
        xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
        yticks(1:numClusters); yticklabels(clusterNames); axis square
        ylabel('Current State'); xlabel('Next New State');
        
        sig_thresh = 0.05 / ((numClusters^2)-numClusters);      % bonferroni correction, for two-tailed p-values so only
        [y,x] = find(reshape(p, 6, 6).*~eye(numClusters) < sig_thresh);
        text(x-.12,y+.12,'*','Color','w');
        
        [~, ~, ~, adj_p] = fdr_bh(reshape(p, 6, 6).*~eye(numClusters));
        [y,x] = find(adj_p < 0.05);
        text(x-.12,y+.12,'x','Color','w');
    
        caxis_bound = max(max(abs(nBackMinusRestTP.*~eye(numClusters))));
        h = colorbar; ylabel(h,[thought_types{2} '-' thought_types{1}]); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
        COLOR_TICK_LABELS(true,true,numClusters);
        title([thought_types{2} '-' thought_types{1}]);
        set(gca,'FontSize',8);
        set(gca,'TickLength',[0 0]);
        set(gca,'Fontname','arial');
    
        if J==0
            sgtitle(["Difference in transition probability between neutral and " thought_type " within control group"])
        elseif J==1
            sgtitle(["Difference in transition probability between neutral and " thought_type " within clinical group"])
        end

        f.PaperUnits = 'inches';
        f.PaperSize = [8 4];
        f.PaperPosition = [0 0 8 4];
        % plot transition probability matrices
        saveas(f,fullfile(savedir,strcat('transitionProbability_', thought_type), '/',['within','diff_', group_type,'_TransProbs_k',num2str(numClusters),'.pdf']),'pdf');
    end
    % see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
    % directed network
    
    %% %% plot transition probabilities group differences (analysis iii)
    % differences between one-back transitions between groups
    % 
    % [~,p1] = ttest2(neutralTransitionProbabilityMatsNoPersist_oneBack(groupInd==0,:),RWTransitionProbabilityNoPersist2D_oneBack(groupInd==1,:));
    % find(p1 < 0.05/numClusters^2)
    % 
    % [~,p2] = ttest2(neutralTransitionProbabilityMatsNoPersist_oneBack(groupInd==0,:),RWTransitionProbabilityNoPersist2D_oneBack(groupInd==0,:));
    % find(p2 < 0.05/numClusters^2)
    % 
    % [~,p3] = ttest2(neutralTransitionProbabilityMatsNoPersist_oneBack(groupInd==1,:),RWTransitionProbabilityNoPersist2D_oneBack(groupInd==1,:));
    % find(p3< 0.05/numClusters^2)
    % 
    % 
    % %diferences from neutral cue to 1-back within groups
    % 
    % [~,p4] = ttest2(neutralTransitionProbabilityMatsNoPersist_thought_cue(groupInd==1,:),neutralTransitionProbabilityNoPersist2D_oneBack(groupInd==1,:));
    % find(p4< 0.05/numClusters^2)
    % 
    % [~,p5] = ttest2(neutralTransitionProbabilityMatsNoPersist_thought_cue(groupInd==0,:),neutralTransitionProbabilityNoPersist2D_oneBack(groupInd==0,:));
    % find(p5< 0.05/numClusters^2)
    % 
    % [~,p6] = ttest2(RWTransitionProbabilityMatsNoPersist_thought_cue(groupInd==1,:),RWTransitionProbabilityNoPersist2D_oneBack(groupInd==1,:));
    % find(p6< 0.05/numClusters^2)
    % 
    % [~,p7] = ttest2(RWTransitionProbabilityMatsNoPersist_thought_cue(groupInd==0,:),RWTransitionProbabilityNoPersist2D_oneBack(groupInd==0,:));
    % find(p7< 0.05/numClusters^2)
    
    
    thought_types = {'Neutral',  char(types(condition+1))};
    group_types = {'control', 'clinical'};
    for J = 0:1
        group_type = group_types{J+1};
        for K = 1:2
            thought_type = thought_types{K};
            
            if strcmpi(thought_type, 'Neutral')
                x = neutralTransitionProbabilityMatsNoPersist_thought_cue(groupInd==J,:);
                y = neutralTransitionProbabilityNoPersist2D_oneBack(groupInd==J,:);
                [~,p] = ttest2(y,x);
                x = reshape(nanmean(x,1), 6, 6);
                y = reshape(nanmean(y,1), 6, 6);
                maxVal = max(max([x,y])); % sync color scales
            else
                x = RWTransitionProbabilityMatsNoPersist_thought_cue(groupInd==J,:);
                y = RWTransitionProbabilityNoPersist2D_oneBack(groupInd==J,:);
                [~,p] = ttest2(y,x);
                x = reshape(nanmean(x,1), 6, 6);
                y = reshape(nanmean(y,1), 6, 6);
                maxVal = max(max([x,y])); % sync color scales
            end
            f = figure;
    
            subplot(1,3,1);
            imagesc(x);
            xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
            xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
            COLOR_TICK_LABELS(true,true,numClusters);
            ylabel('Current State'); xlabel('Next New State');
            title(thought_type);
            set(gca,'FontSize',8);
            set(gca,'TickLength',[0 0]);
            set(gca,'Fontname','arial');
            caxis([0 maxVal]); colorbar
    
            subplot(1,3,2);
            imagesc(y);
            xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
            xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
            COLOR_TICK_LABELS(true,true,numClusters);
    
            ylabel('Current State'); xlabel('Next New State');
            title('1-back');
            set(gca,'FontSize',8);
            set(gca,'TickLength',[0 0]);
            set(gca,'Fontname','arial');
            caxis([0 maxVal]); colorbar
    
            subplot(1,3,3);
            nBackMinusRestTP = (y-x);
            imagesc(nBackMinusRestTP.*~eye(numClusters)); colormap('plasma');
            xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
            yticks(1:numClusters); yticklabels(clusterNames); axis square
            ylabel('Current State'); xlabel('Next New State');
            
            sig_thresh = 0.05 / ((numClusters^2)-numClusters);      % bonferroni correction, for two-tailed p-values so only
            [y,x] = find(reshape(p, 6, 6).*~eye(numClusters) < sig_thresh);
            text(x-.12,y+.12,'*','Color','w');
            
            [~, ~, ~, adj_p] = fdr_bh(reshape(p, 6, 6).*~eye(numClusters));
            [y,x] = find(adj_p < 0.05);
            text(x-.12,y+.12,'x','Color','w');
            
            caxis_bound = max(max(abs(nBackMinusRestTP.*~eye(numClusters))));
            h = colorbar; ylabel(h,['1-back - ' thought_type]); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
            COLOR_TICK_LABELS(true,true,numClusters);
            title(['1-back > ' thought_type]);
            set(gca,'FontSize',8);
            set(gca,'TickLength',[0 0]);
            set(gca,'Fontname','arial');
    
        
            if J==0
                if K==1
                    sgtitle("Transition probability between 1-back and neutral thought within control group")
                elseif K==2 && condition==1
                    sgtitle("Transition probability between 1-back and worry within control group")
                elseif K==2 && condition==2
                    sgtitle("Transition probability between 1-back and rumination within control group")
                end
            elseif J==1
                if K==1
                    sgtitle("Transition probability between 1-back and neutral thought within clinical group")
                elseif K==2 && condition==1
                    sgtitle("Transition probability between 1-back and worry within clinical group")
                elseif K==2 && condition==2
                    sgtitle("Transition probability between 1-back and rumination within clinical group")
                end
            end


            f.PaperUnits = 'inches';
            f.PaperSize = [8 4];
            f.PaperPosition = [0 0 8 4];
            % plot transition probability matrices
            saveas(f,fullfile(savedir,strcat('transitionProbability_', thought_type), '/',['within_', group_type,'_diff_', regexprep(lower(thought_type), ' ', '_'),'ToOneBack_TransProbs_k',num2str(numClusters),'.pdf']),'pdf');
        end
    end
    % see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
    % directed network
    
    %% Plot directed network
    %% load trans probs
    
    % R = load(['results/',name_root,'/analyses/transitionprobabilities/RestCombTransitionProbabilitiesNoPersist_k',...
    %     num2str(numClusters),name_root,'.mat']);
    % restTransProbsFull = reshape(mean(R.transitionProbability,1),[numClusters numClusters])';
    % 
    % N = load(['results/',name_root,'/analyses/nbackblocks/TransProbsNoPersist2back_k',...
    %     num2str(numClusters),name_root,'.mat']);
    % nBackTransProbsFull = reshape(nanmean(N.BlockTransitionProbability,1),[numClusters numClusters])';
    
    % plot
    
    thresh = 0.25;
    % reformat trans prob matrix so that each edge is unidirectional
    restTransProbs = round(grpAvgRest,2,'significant');% - restTransProbsFull';
    restTransProbs(restTransProbs < thresh) = 0;
    R = digraph(restTransProbs);
    
    nBackTransProbs = round(grpAvg2Back,2,'significant');% - nBackTransProbsFull';
    nBackTransProbs(nBackTransProbs < thresh) = 0;
    N = digraph(nBackTransProbs);
    
    all_edges = [R.Edges.Weight;N.Edges.Weight];
    max_edge = max(all_edges(all_edges >0));
    min_edge = min(all_edges(all_edges >0));
    RWidths = 0.1+5*((R.Edges.Weight-min_edge)/max_edge);
    NWidths = 0.1+5*((N.Edges.Weight-min_edge)/max_edge);
    
    f = figure; 
    subplot(1,2,1); h=plot(R,'Layout','circle','LineWidth',RWidths,'EdgeColor',[0 .361 .6335]); axis off; % for text add 'EdgeLabel',R.Edges.Weight,
    labelnode(h,[1 2 3 4 5 6],clusterNames); h.NodeFontSize = 12;
    title('Neutral');
    set(gca,'FontSize',4);
    subplot(1,2,2); h=plot(N,'Layout','circle','LineWidth',NWidths,'EdgeColor',[1 .518 0]); axis off % for text add 'EdgeLabel',N.Edges.Weight,
    labelnode(h,[1 2 3 4 5 6],clusterNames); h.NodeFontSize = 12;
    title('Worry or Rumination');
    set(gca,'FontSize',4);
    f.PaperUnits = 'centimeters';
    f.PaperPosition = [0 0 12 4];
    f.PaperSize = [12 4];
    saveas(f,fullfile(basedir,'results/',['TransProbDigraph_Thresh',num2str(thresh),'_k',num2str(numClusters),'.pdf']));
    
    % plot scale bars for width of lines
    WidthTicks = [0.26 0.3 0.38]; % set gradations in width axis for scale bar
    WidthScale = 0.1+5*((WidthTicks-min_edge)/max_edge);
    f=figure; hold on;
    for j = 1:length(WidthScale)
        line([1 2],[j j],'LineWidth',WidthScale(j),'Color','k')
        text(2.2,j,num2str(WidthTicks(j)),'FontSize',4);
    end
    xlim([0 3]); ylim([0 5]); axis off
    
    f.PaperUnits = 'centimeters';
    f.PaperPosition = [0 0 2 2];
    f.PaperSize = [2 2];
    saveas(f,fullfile(basedir,'results/',['TransProbDigraphWidthAxisScale_Thresh',num2str(thresh),'_k',num2str(numClusters),'.pdf']));
end