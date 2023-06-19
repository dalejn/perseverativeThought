%% Script to apply optimal control
clear
clc

% set directory to main folder with scripts
basedir = '/Users/dalezhou/Desktop/Dropbox/Projects/inProgress/2023-02-perseverantThought/';
cd(basedir)
addpath(genpath('/Users/dalezhou/Desktop/Dropbox/Projects/inProgress/2023-02-perseverantThought/code/'))

% Define numbers and subjects
nreg=214;
nedge=(nreg*(nreg-1))/2;

% Load fMRI task data

yeo_index = dlmread('data/schaefer200x17CommunityAffiliation.txt');
cortical_names = textscan(fopen('data/schaefer200x17NodeNames.txt'), '%s', 'delimiter', '\n');
cortical_names = cortical_names{:};
subcortical_names = textscan(fopen('data/subcortical_order.txt'), '%s', 'delimiter', '\n');
subcortical_names = subcortical_names{:};

% set inputs

savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
distanceMethod = 'correlation'; % distance metric for clustering, we used correlation
numClusters = 6; % set number of clusters -- this requires a separate process of your choosing to select the best number of clusters. See https://arxiv.org/abs/1007.1075 for a useful discussion. 
nreps = 20;	% how many times to repeat clustering. will choose lowest error solution

nNodes = size(nreg, 1);

% get consistent order of regions from structural connectivity matrix

% Read in structural connectivity matrices
sc_network_files = dir('data/struct_MURI/*.mat');
nfiles = length(sc_network_files);
fa_sq = zeros(nfiles, nedge);
subj_name_vector = {};

for k=1:nfiles
    load(['data/struct_MURI/' sc_network_files(k).name])
    struc_subcortical_order = cellstr(schaefer214x17_region_labels);
    [~, idx_reorder_subcortical_nodes] = ismember(subcortical_names, {struc_subcortical_order{201:214}}');
    idx_reorder_subcortical_nodes = nonzeros(idx_reorder_subcortical_nodes);
    schaefer214x17_count_pass_connectivity(201:214);
    schaefer214x17_count_pass_connectivity2 = schaefer214x17_count_pass_connectivity([1:200 idx_reorder_subcortical_nodes'+200],[1:200 idx_reorder_subcortical_nodes'+200]);
    fa_sq(k,:) = squareform(schaefer214x17_count_pass_connectivity2);
    subj_name = split(sc_network_files(k).name, '-');
    subj_name = split(subj_name{2}, '_');
    subj_name = lower(subj_name{1})
    subj_name_vector{k, 1} = subj_name;
end

subj_name_vector = strrep(subj_name_vector, 'murip','muri')

average_structural_network = squareform(mean(fa_sq,1));
% dlmwrite('data/struct_MURI/average_structural_network.txt', average_structural_network, ' ');

% the order of subcortical nodes in structural network
order_struct = schaefer214x17_region_labels(idx_reorder_subcortical_nodes'+200,:);

% the order of subcortical nodes in time series
order_func = subcortical_names(find(ismember(subcortical_names, {struc_subcortical_order{201:214}})));

%% load time series and get BOLD data formatted as a T-by-nparc matrix
% gather subject dir
subj_ID = dir('data/timeSeries/*_sub-*');
missingFiles = {};

types = {'neutral', 'worry', 'rumination'}

nTR_all = [];
TS = [];
scanInd = [];
nbackInd = [];
groupInd = [];

for k = 1:length(subj_ID)
    currentSub = subj_ID(k).name
    nTR = 0;
    group_type = split(currentSub, '_');
    group_type = group_type(1);
    if strcmpi(group_type, 'control')
        groupInd = [groupInd; 0];
    else
        groupInd = [groupInd; 1];
    end
    for typenum = 1:length(types)
        type = types{typenum};
        % Read in task fMRI
        task_files = dir(['data/timeSeries/' currentSub '/' type '/' type(1) '_ts*.txt']);

        for j = 1:length(task_files)
            try
                task_ts = dlmread(['data/timeSeries/' currentSub '/' type '/' task_files(j).name]);
                task_ts = task_ts([1:200 find(ismember(subcortical_names, {struc_subcortical_order{201:214}}))'+200],:);
                task_ts = task_ts';
                nTR_current = size(task_ts,1);
                nTR = nTR + nTR_current;
                if strcmpi(type, 'neutral')
                    scanInd_current = false(nTR_current,1);
                else
                    scanInd_current = true(nTR_current,1);
                end
                nbackInd_current = [repelem(0,10) repelem(1,nTR_current-10)]';
                TS = [TS; task_ts];
                scanInd = [scanInd; scanInd_current];
                nbackInd = [nbackInd; nbackInd_current];
                

            catch missingFile
                missingFiles = [missingFiles; ['data/timeSeries/' currentSub '/' type '/' task_files(j).name]];
            end
        end
    end
    nTR_all = [nTR_all; nTR];
end

%% load BOLD data

% replace TS with your BOLD data formatted as a T-by-nparc matrix
% where T is the number of time points and nparc is the number of ROIs
D = load('data/concTS.mat'); % HCP U100 sample in lausanne250 parcellation
TS = D.TS;

nTR = sum(nTR_all); % use sa
nsubjs = 44;		% number of subjects is total # of TRs divided by # of TRs per subjects
subjInd = [repelem(1:nsubjs,nTR_all)]'; % index data from each subject
% scanInd = scanInd; % index scans within your time series data
scanlab = {'neutral','worry_or_ruminate'}; % text labels for scan indices
[T,nparc] = size(TS);

%% identify brain states
[partition,~,sumd] = kmeans(TS,numClusters,'Distance', distanceMethod,'Replicates',nreps);

%save('results/partition.mat', 'partition')

% see supplement for information on evaluating choices for the number of clusters
% in our experience, 4-6 clusters are around the elbow of the variance
% explained curve and are robust to subsampling.