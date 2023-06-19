%% control energy

A= dlmread('data/struct_MURI/average_structural_network.txt'); % load example group average structural A matrix -- this is from PNC while fMRI is HCP so DTI is younger than fMRI here
c = 0; T = 5; % set time scale parameters based on values from paper
Anorm = NORMALIZE(A,c); % normalize A by maximum eigenvalue - eye(N) to make marginally stable

% define x0 and xf, initial and final states as cluster centroids for each
% state transition
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities

x0 = centroids(:,Xo_ind);
xf = centroids(:,Xf_ind); % now each column of x0 and xf represent state transitions
WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
E_full = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition

% compute weighted control energy:
network7labels = dlmread(fullfile(basedir,'data','schaefer200x17CommunityAffiliation_binned.txt'));
InputVector = ismember(network7labels(1:nparc),1); % weight input towards visual system as in task
B = InputVector .*eye(nparc) + eye(nparc); % construct input matrix allowing input only into selected regions in InputVector
E_weighted = zeros(1,numClusters^2);
for transition = 1:numClusters^2    
    [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
    E_weighted(transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
end
