%% plot brain states
hexMap = GET_CLUSTER_COLORS(6);
for k = 1 : length(hexMap)
	thisCell = hexMap{k}
	r = hex2dec(thisCell(1:2))
	g = hex2dec(thisCell(3:4))
	b = hex2dec(thisCell(5:6))
	myColorMap(k, :) = [r, g, b]
end
myColorMap = myColorMap / 255; % Normalize to range 0-1

%find row of last clinical subj
row_ind = find(groupInd);

fig=figure(1);
hold on
for N = 1:nsubjs
    subplot(nsubjs,1,N);
	subjPartition = partition(subjInd == N);
    
%     % to figure out where to draw group boundary line
%     if N==row_ind(end)
%         imagesc(repelem(1, length(subjPartition)))
%         continue
%     end
    imagesc(subjPartition')
    colormap(myColorMap)
    set(gca,'XTick',[]); %which will get rid of all the markings for the y axis
    set(gca,'YTick',[]); %which will get rid of all the markings for the y axis
    set(gca,'Yticklabel',[]) 
    set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
end
hold off
h = axes(fig,'visible','off'); 
c = lcolorbar(["1" "2" "3" "4" "5" "6"]);
set(c, 'Position', [0.93 0.168 0.022 0.7])
annotation('line', [0.115 0.925], [0.4075 0.4075], 'Color', 'r', 'LineWidth', 1);
annotation('textbox',[0 0.1 0.6 0.6],'String', {'Clinical group'}, 'EdgeColor','none')
annotation('textbox',[0 0.1 0.2 0.2],'String', {'Control group'}, 'EdgeColor','none')

%% compute centroids and plot
centroids = GET_CENTROIDS(TS,partition,numClusters);
% name clusters based on alignment with Yeo resting state networks
clusterNames = NAME_CLUSTERS_ANGLE(centroids);  % need to add a prior Yeo partition labels for your parcellation
[clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(centroids);  % need to add a prior Yeo partition labels for your parcellation

f = figure;
subplot(1,2,1); imagesc(centroids); title('Centroids'); xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
colormap('plasma'); axis square; colorbar; set(gca,'FontSize',8); COLOR_TICK_LABELS(true,false,numClusters);
subplot(1,2,2); imagesc(corr(centroids)); title('Centroid Similarity'); colorbar; caxis([-1 1]); 
colormap('plasma'); axis square; set(gca,'FontSize',8); xticks(1:numClusters); yticks(1:numClusters); 
xticklabels(clusterNames); yticklabels(clusterNames); xtickangle(90);
COLOR_TICK_LABELS(true,true,numClusters);
f.PaperUnits = 'inches';
f.PaperSize = [4 2];
f.PaperPosition = [0 0 4 2];
saveas(f,fullfile(savedir,['Centroids_k',num2str(numClusters),'.pdf']));

% can use python scripts to make surface plots