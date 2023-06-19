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

%% generate all permutations to reorder columns of k-means clustering partition

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
