clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Predict state swiching for tTIS patients based on the WEiDA State Space
% model established in DBS

%%%%%%

addpath('./Functions')

Vemp = csvread("./WEiDA4_atlasAAL78_DBS/Vemp.csv");    % Cluster centroid weighted-eigenvector (PC pattern) for each brain state in DBS
Vemp = Vemp';
NumClusters = 4;    % Number of brain state in DBS state space
TR = 2.25;       % Repetition Time (seconds)
N_ti = 11;      % number of tTIS subjects

DirOut='./WEiDA4_atlasAAL78_tTIS/';
mkdir(DirOut);

pth = './data/tTIS/';

for r=[1,3]      % Condition of 1 for tTIS-PRE or 3 for tTIS-ACT
    rawfile_run = dir(['./data/tTIS/sub-tTIS*_run-0' num2str(r) '_BOLD_AAL78_MNI_basic36_sm6.csv']);
    for nsub=1:N_ti
        file = rawfile_run(nsub).name;
        subpth = [pth,file];
        a = importdata(subpth);
        a = a.data;
        a(:,1)=[];
        BOLD = a';      % Bold signal for a subject under a condition of tTIS-PRE or tTIS-POST (N_areas,TR)
        % Predict state switching based on DBS state space
        % State fractional occupancy, State dwell time, transition
        % probability, State time series
        [PTR,Pstates(nsub,:),LTime(nsub,:),TsState(nsub,:)]=WEiDA_fix_cluster(BOLD,NumClusters,Vemp,TR);  
        csvwrite([DirOut 'Ptrans' num2str(11*(r-1)+nsub) '.csv'],PTR);
    end
    csvwrite([DirOut 'LT' num2str(r) 'emp.csv'],LTime);
    csvwrite([DirOut 'P' num2str(r) 'emp.csv'],Pstates);
    csvwrite([DirOut 'State_ts' num2str(r) '.csv'],TsState);
end
