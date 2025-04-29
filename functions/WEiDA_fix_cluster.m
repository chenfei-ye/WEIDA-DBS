function [PTRANSITION,Pstates,LTime,IDX] = WEiDA_fix_cluster(BOLDsig,NumClusters,Center,TR)

[N_areas, Tmax]=size(BOLDsig);

% Preallocate variables to save PC patterns and associated information
Weighted_Eig=zeros(Tmax,1*N_areas); % All weighted eigenvectors

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)

Phase_BOLD=zeros(N_areas,Tmax);

% Get the BOLD phase using the Hilbert transform
for seed=1:N_areas
    BOLDsig(seed,:)=demean(detrend(BOLDsig(seed,:)));
    signal_filt =filtfilt(bfilt,afilt,BOLDsig(seed,:));
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
    Weighted_Eig(t_all,:)=V1; %vertcat(V1,V2);
end

clear signal_filt iPC V1 Phase_BOLD

%% 2 - Predict state for each TR (Weighted eigenvector)
IDX=zeros(t_all,1);

for t=1:t_all
    for j=1:NumClusters
        di(j)=sqrt(sum((Weighted_Eig(t,:)-Center(j,:)).^2));  % Euclidean distance between Weighted eigenvector and State PC patterns
    end
    [aux indmin]=min(di);   % Determine the brain state for each TR (Weighted eigenvector) based on Euclidean distance
    IDX(t)=indmin;
end

Pstates=zeros(1,NumClusters);
for c=1:NumClusters
    Pstates(c)=mean(IDX==c);
    % Mean Lifetime
    Ctime_bin=IDX==c;
    
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
        C_Durations=0;
    end
    % State dwell time
    LTime(c)=mean(C_Durations)*TR;
end

% State fractional occupancy
Pstates=Pstates/sum(Pstates);

% Transition probability
PTRANSITION=zeros(NumClusters,NumClusters);
i=1;
for c1=1:NumClusters
    j=1;
    for c2=1:NumClusters
        sumatr=0;
        for t=1:length(IDX)-1
            if IDX(t)==c1 && IDX(t+1)==c2
                sumatr=sumatr+1;
            end
        end
        if length(find(IDX(1:length(IDX)-1)==c1)) ~= 0
            PTRANSITION(i,j)=sumatr/(length(find(IDX(1:length(IDX)-1)==c1)));
        end
        j=j+1;
    end
    i=i+1;
end

 