%%
load('results/partition.mat', 'partition')

%% test clustering on splits

rng(1);
nreps = 500;
numClusters = 6;

concTS = TS;
%nobs = size(TS,1);
nobs = max(subjInd);

disp('start k-means');
for R = 1:nreps
	obs1 = randperm(nobs,floor(0.5*nobs));
    [partition1,~,sumd1] = kmeans(concTS(obs1,:),numClusters,'Distance',distanceMethod, 'MaxIter',500);
    partition1 = int8(partition1);

    obs2 = find(~ismember(1:nobs,obs1));	% cluster the other half
    [partition2,~,sumd2] = kmeans(concTS(obs2,:),numClusters,'Distance',distanceMethod, 'MaxIter',500);
    partition2 = int8(partition2); 
    save(['SubjSHkmeans',num2str(R),'k_',num2str(numClusters),'rep',num2str(R),'.mat'],'partition1','sumd1','obs1','partition2','sumd2','obs2');
    clear partition1; clear partition2;
    clear sumd1; clear sumd2;
    disp(['K-means ',num2str(R)]);
end

disp('complete');

%%

r_hist = [];
splits_files = dir('SubjSHkmeans*6rep*.mat');
for file_index=1:nreps
%for file_index=1:1
    load(splits_files(file_index).name)
    centroids_partition1 = GET_CENTROIDS(TS,partition1,numClusters);
    centroids_partition2 = GET_CENTROIDS(TS,partition2,numClusters);
%     findMax = corr(centroids_partition1, centroids_partition2);
%     reorderedIndex = [];
%     for k=1:numClusters
%        test_value = 2;
%        [M,I] = max(findMax(k,:));       
%        % for ties
%        while ismember(I, reorderedIndex)
%          [mx,ind]=maxk(findMax(k,:),test_value);
%          I = ind(test_value);
%          test_value = test_value + 1;
%        end
%        
%        reorderedIndex = [reorderedIndex I]; 
%     end
    A = centroids_partition1;
    B = centroids_partition2;
    numColumns = size(A, 2);
    permutations = perms(1:numColumns);  % Generate all permutations of column indices

    maxCorrelation = 0;
    bestPermutation = [];

    for i = 1:size(permutations, 1)
        permutedB = B(:, permutations(i, :));  % Permute columns of A
        correlation = corrcoef(A, permutedB);  % Calculate correlation with B

        if correlation(2,1) > maxCorrelation
            maxCorrelation = correlation(2,1);
            bestPermutation = permutations(i, :);
        end
    end

    reorderedB = B(:, bestPermutation);  % Reorder columns of A based on best permutation
    
    [r,p] = corrcoef(A, reorderedB);
    r_hist = [r_hist r(2)];
end
hist(r_hist)

%% split by 1-back or thought

rng(1);
nreps = 1;
numClusters = 6;

concTS = TS;
%nobs = size(TS,1);
nobs = max(subjInd);

disp('start k-means');
for R = 1:nreps;
    [partition1,~,sumd1] = kmeans(concTS(nbackInd==0,:),numClusters,'Distance',distanceMethod, 'MaxIter',500);
    partition1 = int8(partition1);

    [partition2,~,sumd2] = kmeans(concTS(nbackInd==1,:),numClusters,'Distance',distanceMethod, 'MaxIter',500);
    partition2 = int8(partition2); 
end

centroids_partition1 = GET_CENTROIDS(TS,partition1,numClusters);
centroids_partition2 = GET_CENTROIDS(TS,partition2,numClusters);
% findMax = corr(centroids_partition1, centroids_partition2);
% reorderedIndex = [];
% for k=1:numClusters
%    test_value = 2;
%    [M,I] = max(findMax(k,:));       
%    % for ties
%    while ismember(I, reorderedIndex)
%      [mx,ind]=maxk(findMax(k,:),test_value);
%      I = ind(test_value);
%      test_value = test_value + 1;
%    end
%    reorderedIndex = [reorderedIndex I]; 
% end
% [r,p] = corrcoef(centroids_partition1, centroids_partition2(:,reorderedIndex));
% 
% similarityMatrix = corr(centroids_partition1, centroids_partition2(:,reorderedIndex));
A = centroids_partition1;
B = centroids_partition2;
numColumns = size(A, 2);
permutations = perms(1:numColumns);  % Generate all permutations of column indices

maxCorrelation = 0;
bestPermutation = [];

for i = 1:size(permutations, 1)
    permutedB = B(:, permutations(i, :));  % Permute columns of A
    correlation = trace(corr(A, permutedB));  % Calculate correlation with B

    if correlation > maxCorrelation
        maxCorrelation = correlation;
        bestPermutation = permutations(i, :);
    end
end

reorderedB = B(:, bestPermutation);  % Reorder columns of A based on best permutation

similarityMatrix = corr(A, reorderedB);

heatmap(similarityMatrix);
colormap(jet(512))

%% repeat clustering and find centroid name distribution for each state

centroids_partition1 = GET_CENTROIDS(TS,partition,numClusters);
clusterNames = NAME_CLUSTERS_ANGLE(centroids_partition1);  % need to add a prior Yeo partition labels for your parcellation

names_all = cell(nreps,numClusters);

% identify brain states
for R = 1:500
    [partition_new,~,sumd_new] = kmeans(TS,numClusters,'MaxIter',nreps, 'Distance', distanceMethod,'Replicates',5);
    
    centroids_partition2 = GET_CENTROIDS(TS,partition_new,numClusters);
    
    A = centroids_partition1;
    B = centroids_partition2;
    numColumns = size(A, 2);
    permutations = perms(1:numColumns);  % Generate all permutations of column indices

    maxCorrelation = 0;
    bestPermutation = [];

    for i = 1:size(permutations, 1)
        permutedB = B(:, permutations(i, :));  % Permute columns of A
        correlation = trace(corr(A, permutedB));  % Calculate correlation with B

        if correlation > maxCorrelation
            maxCorrelation = correlation;
            bestPermutation = permutations(i, :);
        end
    end

    reorderedB = B(:, bestPermutation);  % Reorder columns of B based on best permutation
    
    clusterNames_new = NAME_CLUSTERS_ANGLE(reorderedB);  % need to add a prior Yeo partition labels for your parcellation
    names_all(R,:) = clusterNames_new;

    save(['repClustering/partition_' num2str(R) '.mat'], 'partition_new', 'clusterNames_new')
end

%% plot pie charts

t = tiledlayout(1,6);
nexttile
pie(categorical(names_all(:,1)))
nexttile
pie(categorical(names_all(:,2)))
nexttile
pie(categorical(names_all(:,3)))
nexttile
pie(categorical(names_all(:,4)))
nexttile
pie(categorical(names_all(:,5)))
nexttile
pie(categorical(names_all(:,6)))

%% plot spider charts

% VIS SOM DAT VAT LIM FPN DMN
% 1 94.2 5.8 0 0 0 0 0
% 2 69.8 25 0 4.4 0 0.8 0
% 3 21.4 12.6 0 0 0 66 0
% 4 23.2 76.8 0 0 0 0 0
% 5 0 20.6 0 0 0 0 79.4
% 6 0 84 4.8 0.2 0 0 11

% Initialize data points
D1 = [94.2 5.8 0 0 0 0 0];
D2 = [69.8 25 0 4.4 0 0.8 0];
D3 = [21.4 12.6 0 0 0 66 0];
D4 = [23.2 76.8 0 0 0 0 0];
D5 = [0 20.6 0 0 0 0 79.4];
D6 = [0 84 4.8 0.2 0 0 11];

P = [D1; D2; D3; D4; D5; D6];

% Spider plot
spider_plot(P,...
    'AxesLabels', {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'},...
    'AxesLimits', [0, 0, 0, 0, 0, 0, 0;100, 100, 100, 100, 100, 100, 100],... % [min axes limits; max axes limits]
    'AxesInterval', 2,...
    'FillOption', {'on', 'on', 'on', 'on', 'on', 'on'},...
    'FillTransparency', [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]);


%% function to generate all permutations to reorder columns of k-means clustering partition

function reorderedA = maximizeCorrelation(A, B)
    numColumns = size(A, 2);
    permutations = perms(1:numColumns);  % Generate all permutations of column indices
    
    maxCorrelation = -Inf;
    bestPermutation = [];
    
    for i = 1:size(permutations, 1)
        permutedA = A(:, permutations(i, :));  % Permute columns of A
        correlation = corr(permutedA, B);  % Calculate correlation with B
        
        if correlation > maxCorrelation
            maxCorrelation = correlation;
            bestPermutation = permutations(i, :);
        end
    end
    
    reorderedA = A(:, bestPermutation);  % Reorder columns of A based on best permutation
end