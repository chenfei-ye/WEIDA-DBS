clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate state metastability for each subject

%%%%%%

addpath('./functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data preppare
pth_stn = './data/DBS/';
pth_ti = './data/tTIS/';
rawfile_hc = dir("./data/DBS/sub-HC*_BOLD_AAL78_MNI_basic36_sm6.csv");
rawfile_stn = dir("./data/DBS/sub-STN*_BOLD_AAL1+PD25_MNI_basic36_sm6.csv");
rawfile_tipre = dir("./data/tTIS/sub-*_active_run-01_BOLD_AAL78_MNI_basic36_sm6.csv");
rawfile_tipost = dir("./data/tTIS/sub-*_active_run-03_BOLD_AAL78_MNI_basic36_sm6.csv");
state_stn = csvread("./WEiDA4_atlasAAL78_DBS/State_ts.csv");
state_ti1 = csvread("./WEiDA4_atlasAAL78_tTIS/State_ts1.csv");
state_ti3 = csvread("./WEiDA4_atlasAAL78_tTIS/State_ts3.csv");
state_ti = [state_ti1;state_ti3];

Vemp = csvread("./WEiDA4_atlasAAL78_DBS/Vemp.csv");

N_hc = 16;
N_stn = 25;
N_ti = 11;

for i=1:N_hc
    file = rawfile_hc(i).name;
    subpth = [pth_stn,file];
    a = importdata(subpth);
    a = a.data;
    a(:,1)=[];
    ts_cell{i,1} = a';
end
    
for i=1:2*N_stn
    file = rawfile_stn(i).name;
    subpth = [pth_stn,file];
    a = importdata(subpth);
    a = a.data;
    a(:,1)=[];
    ts_cell{N_hc+i,1} = a';
end

for i=1:N_ti
    file = rawfile_tipre(i).name;
    subpth = [pth_ti,file];
    a = importdata(subpth);
    a = a.data;
    a(:,1)=[];
    ts_cell1{i,1} = a';
end

for i=1:N_ti
    file = rawfile_tipost(i).name;
    subpth = [pth_ti,file];
    a = importdata(subpth);
    a = a.data;
    a(:,1)=[];
    ts_cell1{N_ti+i,1} = a';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate state metastability for HC and DBS

[n_Subjects, ~]=size(ts_cell);
[N_areas, Tmax]=size(ts_cell{1,1});
Weighted_Eig=zeros(Tmax,1*N_areas); % All weighted eigenvectors

% Parameters of the data
TR=2.1;  % Repetition Time (seconds)

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
    % Get the BOLD signals from this subject in this condition
    BOLD = ts_cell{s,1};
    Phase_BOLD=zeros(N_areas,Tmax);
    state_ts=state_stn((s-1)*Tmax+1:s*Tmax);

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
        Weighted_Eig(t,:)=V1; %vertcat(V1,V2);
    end
    
    % Calculate VAR metastability for each state
    for st=1:4
        OPs(s,st)=mean(sum((Weighted_Eig(state_ts==st,:)-Vemp(:,st)').^2) / (sum(state_ts==st)-1));
    end
end

csvwrite("./WEiDA4_atlasAAL78_DBS/Metastability_HC_DBS.csv",OPs);

%% Calculate state metastability for tTIS

[n_Subjects, ~]=size(ts_cell1);
[N_areas, Tmax]=size(ts_cell1{1,1});
Weighted_Eig=zeros(Tmax,1*N_areas); % All weighted eigenvectors

% Parameters of the data
TR=2.25;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

for s=1:n_Subjects
    [N_areas, Tmax]=size(ts_cell1{s,1});
    % Get the BOLD signals from this subject in this condition
    BOLD = ts_cell1{s,1};
    Phase_BOLD=zeros(N_areas,Tmax);
    state_ts=state_ti(s,:);

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
        Weighted_Eig(t,:)=V1; %vertcat(V1,V2);
    end
    
    % Calculate VAR metastability for each state
    for st=1:4
        OPs1(s,st)=mean(sum((Weighted_Eig(state_ts==st,:)-Vemp(:,st)').^2) / (sum(state_ts==st)-1));
    end
end

csvwrite("./WEiDA4_atlasAAL78_tTIS/Metastability_tTIS.csv",OPs1);

    