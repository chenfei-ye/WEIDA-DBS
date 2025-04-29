clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WEiDA state space modelling pipeline
% (Based on DBS data)

%%%%%%

addpath('./functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data preppare

pth = './data/DBS/';
rawfile_hc = dir("./data/DBS/sub-HC*_BOLD_AAL78_MNI_basic36_sm6.csv");
rawfile_stn = dir("./data/DBS/sub-STN*_BOLD_AAL78_MNI_basic36_sm6.csv");

N_hc = 16;
N_stn = 25;

for i=1:N_hc
    file = rawfile_hc(i).name;
    subpth = [pth,file];
    a = importdata(subpth);
    a = a.data;
    a(:,1)=[];
    ts_cell{i,1} = a';
    ts_cell{i,2} = 1;
end
    

for i=1:2*N_stn
    file = rawfile_stn(i).name;
    subpth = [pth,file];
    a = importdata(subpth);
    a = a.data;
    a(:,1)=[];
    ts_cell{N_hc+i,1} = a';
end
for i=1:N_stn
    ts_cell{N_hc+i,2} = 2;    
end
for i=N_stn+1:2*N_stn
    ts_cell{N_hc+i,2} = 3;    
end
n_cond=3;   % Number of conditions: 1 for HC, 2 for DBS-OFF, 3 for DBS-ON

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 - Compute the Weighted Eigenvectors from the BOLD datasets
disp('Processing the eigenvectors from BOLD data')
% Load here the BOLD data (which may be in different formats)
% Here the BOLD time courses in AAL parcellation are organized as cells,
% where ts_cell{1,1} corresponds to the BOLD data from subject 1 in
% condition 1 and contains a matrix with lines=N_areas and columns=Tmax.

[n_Subjects, ~]=size(ts_cell);
[N_areas, Tmax]=size(ts_cell{1,1});

% Parameters of the data
TR=2.1;  % Repetition Time (seconds)

% Preallocate variables to save PC patterns and associated information
WEiDA_Eig=zeros(Tmax*n_Subjects,1*N_areas); % All Weighted eigenvectors
Time_all=zeros(2, n_Subjects*Tmax); % vector with subject nr and cond at each t
t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

for s=1:n_Subjects
    [N_areas, Tmax]=size(ts_cell{s,1});
    cond = ts_cell{s,2};
    % Get the BOLD signals from this subject in this condition
    BOLD = ts_cell{s,1};
    Phase_BOLD=zeros(N_areas,Tmax);

    % Get the BOLD phase using the Hilbert transform
    for seed=1:N_areas
        ts=demean(detrend(BOLD(seed,:)));
        signal_filt =filtfilt(bfilt,afilt,ts);
        Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
    end

    for t=1:Tmax

        %Calculate the Instantaneous PC (BOLD Phase Synchrony)
        iPC=zeros(N_areas);
        for n=1:N_areas
            for p=1:N_areas
                iPC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
            end
        end

        % Get the Weighted eigenvector

        [V,lamda]=eigs(iPC);
        V1=zeros(N_areas,1);
        for i=1:size(V,2)
            V1=V1+lamda(i,i)*V(:,i);
        end
        V1=V1/sum(sum(lamda));
        % Make sure the largest component is negative
        if mean(V1>0)>.5
            V1=-V1;
        elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
            V1=-V1;
        end           
        % Save V1 from all frames in all fMRI sessions in Weighted eig
        t_all=t_all+1; % Update time
        WEiDA_Eig(t_all,:)=V1; %vertcat(V1,V2);
        Time_all(:,t_all)=[s cond]; % Information that at t_all, V1 corresponds to subject s in a given condition   
    end
end
clear BOLD tc_aal signal_filt iPC VV V1 V2 Phase_BOLD

%% 2 - Cluster the WEiDA Eigenvectors

disp('Clustering the eigenvectors into')
% Weighted_Eig is a matrix containing all the eigenvectors:
% Collumns: N_areas are brain areas (variables)
% Rows: Tmax*n_Subjects are all time points (independent observations)

% Set maximum/minimum number of clusters
% There is no fixed number of states the brain can display
% Extending depending on the hypothesis of each work
mink=4;
maxk=4;
rangeK=mink:maxk;

% Set the parameters for Kmeans clustering
Kmeans_results=cell(size(rangeK));
for k=1:length(rangeK)
    disp(['- ' num2str(rangeK(k)) ' clusters'])
    [IDX, C, SUMD, D]=kmeans(WEiDA_Eig,rangeK(k),'Distance','sqeuclidean','Replicates',300,'Display','off','MaxIter',500);
    Kmeans_results{k}.IDX=IDX;   % Cluster time course - numeric collumn vectos
    Kmeans_results{k}.C=C;       % Cluster centroids (PC patterns)
    Kmeans_results{k}.SUMD=SUMD; % Within-cluster sums of point-to-centroid distances
    Kmeans_results{k}.D=D;       % Distance from each point to every centroid
    ss=silhouette(WEiDA_Eig,IDX,'sqeuclidean');
    sil(k)=mean(ss);
end

%% 3 - Analyze clustering results

% For every fMRI scan calculate probability and lifetimes of each pattern c.
P=zeros(n_cond,n_Subjects,maxk-mink+1,maxk);
LT=zeros(n_cond,n_Subjects,maxk-mink+1,maxk);

for k=1:length(rangeK)
    for cond=1:n_cond   
        for s=1:n_Subjects
            
            % Select the time points representing this subject and task
            T=((Time_all(1,:)==s)+(Time_all(2,:)==cond))>1;
            Ctime=Kmeans_results{k}.IDX(T);
            
            for c=1:rangeK(k)
                % Probability
                P(cond,s,k,c)=mean(Ctime==c);
                
                % Mean Lifetime
                Ctime_bin=Ctime==c;
                
                % Detect switches in and out of this state
                a=find(diff(Ctime_bin)==1);
                b=find(diff(Ctime_bin)==-1);
                
                % We discard the cases where state sarts or ends ON
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end
                if ~isempty(a) && ~isempty(b)
                    C_Durations=b-a;
                else
                    C_Durations=NaN;
                end
                LT(cond,s,k,c)=mean(C_Durations)*TR;
            end
        end
    end
end

%% 4 - Calculate state dynamics

disp(' ')
disp('%%% PLOTS %%%%')
disp(['Choose number of clusters between ' num2str(rangeK(1)) ' and ' num2str(rangeK(end)) ])

K = input('Number of clusters: ');
Number_Clusters=K;
Best_Clusters=Kmeans_results{rangeK==K};
k=find(rangeK==K);

% Clusters are sorted according to their probability of occurrence
ProbC=zeros(1,K);
for c=1:K
    ProbC(c)=mean(Best_Clusters.IDX==c);
end
[~, ind_sort]=sort(ProbC,'descend');

% Get the K patterns
V=Best_Clusters.C(ind_sort,:);
[~, N]=size(Best_Clusters.C);

Vemp=V;   % Cluster centroid weighted-eigenvector (PC pattern) for each brain state

% State fractional occupancy
P1emp=squeeze(P(1,:,k,ind_sort));
P2emp=squeeze(P(2,:,k,ind_sort));
P3emp=squeeze(P(3,:,k,ind_sort));
% State dwell time
LT1emp=squeeze(LT(1,:,k,ind_sort));
LT2emp=squeeze(LT(2,:,k,ind_sort));
LT3emp=squeeze(LT(3,:,k,ind_sort));

PTR1emp=zeros(K);
PTR2emp=zeros(K);
PTR3emp=zeros(K);
n_sub1=zeros(K,1);
n_sub2=zeros(K,1);
n_sub3=zeros(K,1);

% Average transition probability for each condition (1 for HC, 2 for DBS-OFF, 3 for DBS-ON)
for cond=1:3 
    for s=1:n_Subjects
        % Select the time points representing this subject and cond
        T=((Time_all(1,:)==s)+(Time_all(2,:)==cond))>1;
        Ctime=Kmeans_results{k}.IDX(T);

        if cond==1
            i=1;
            for c1=ind_sort
                j=1;
                for c2=ind_sort
                    sumatr=0;
                    for t=1:length(Ctime)-1
                        if Ctime(t)==c1 && Ctime(t+1)==c2
                            sumatr=sumatr+1;
                        end
                    end
                    if length(find(Ctime(1:length(Ctime)-1)==c1)) ~= 0
                        PTR1emp(i,j)=PTR1emp(i,j)+sumatr/(length(find(Ctime(1:length(Ctime)-1)==c1)));
                    end
                    j=j+1;
                end
                if length(find(Ctime(1:length(Ctime)-1)==c1)) ~=0
                    n_sub1(i)=n_sub1(i)+1;
                end
                i=i+1;
            end
        end

        if cond==2
            i=1;
            for c1=ind_sort
                j=1;
                for c2=ind_sort
                    sumatr=0;
                    for t=1:length(Ctime)-1
                        if Ctime(t)==c1 && Ctime(t+1)==c2
                            sumatr=sumatr+1;
                        end
                    end
                    if length(find(Ctime(1:length(Ctime)-1)==c1)) ~=0
                        PTR2emp(i,j)=PTR2emp(i,j)+sumatr/(length(find(Ctime(1:length(Ctime)-1)==c1)));
                    end
                    j=j+1;
                end
                if length(find(Ctime(1:length(Ctime)-1)==c1)) ~=0
                    n_sub2(i)=n_sub2(i)+1;
                end

                i=i+1;
            end
        end

        if cond==3
            i=1;
            for c1=ind_sort
                j=1;
                for c2=ind_sort
                    sumatr=0;
                    for t=1:length(Ctime)-1
                        if Ctime(t)==c1 && Ctime(t+1)==c2
                            sumatr=sumatr+1;
                        end
                    end
                    if length(find(Ctime(1:length(Ctime)-1)==c1)) ~=0
                        PTR3emp(i,j)=PTR3emp(i,j)+sumatr/(length(find(Ctime(1:length(Ctime)-1)==c1)));
                    end
                    j=j+1;
                end
                if length(find(Ctime(1:length(Ctime)-1)==c1)) ~=0
                    n_sub3(i)=n_sub3(i)+1;
                end

                i=i+1;
            end
        end

    end
end
for i=1:K
    PTR1emp(i,:)=PTR1emp(i,:)/n_sub1(i);
    PTR2emp(i,:)=PTR2emp(i,:)/n_sub2(i);
    PTR3emp(i,:)=PTR3emp(i,:)/n_sub3(i);   
end

% Calculate transition probability for each subject
DirOut=['./WEiDA' num2str(K) '_atlasAAL78_DBS/'];
mkdir(DirOut);

for s=1:n_Subjects
    PTRemp=zeros(K);
    % Select the time points representing this subject
    T=(Time_all(1,:)==s);
    Ctime=Best_Clusters.IDX(T);
    i=1;
    for c1=ind_sort
        j=1;
        for c2=ind_sort
            sumatr=0;
            for t=1:length(Ctime)-1
                if Ctime(t)==c1 && Ctime(t+1)==c2
                    sumatr=sumatr+1;
                end
            end
            if length(find(Ctime(1:length(Ctime)-1)==c1)) ~= 0
                PTRemp(i,j)=sumatr/(length(find(Ctime(1:length(Ctime)-1)==c1)));
            end
            j=j+1;
        end
        i=i+1;
    end
    csvwrite([DirOut 'Ptrans' num2str(s) '.csv'],PTRemp);
end

csvwrite([DirOut 'Vemp.csv'],Vemp.');
csvwrite([DirOut 'LT1emp.csv'],LT1emp);
csvwrite([DirOut 'LT2emp.csv'],LT2emp);
csvwrite([DirOut 'LT3emp.csv'],LT3emp);
csvwrite([DirOut 'P1emp.csv'],P1emp);
csvwrite([DirOut 'P2emp.csv'],P2emp);
csvwrite([DirOut 'P3emp.csv'],P3emp);
csvwrite([DirOut 'PTR1emp.csv'],PTR1emp);
csvwrite([DirOut 'PTR2emp.csv'],PTR2emp);
csvwrite([DirOut 'PTR3emp.csv'],PTR3emp);

% State switching (State time series)
IDX = Kmeans_results{k}.IDX;
for i=1:K
    IDX(Best_Clusters.IDX==i)=find(ind_sort==i);
end

csvwrite([DirOut 'State_ts.csv'],IDX);



