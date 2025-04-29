clc; clear all;
% Calculate overlap with the 7 resting-state Networks from Yeo et al. 2011

addpath('./functions')
DirOut = './WEiDA4_atlasAAL78_DBS/';

% load the masks (AAL and Yeo networks) in MNI space 2mm
file1 = gunzip("./data/atlas/AAL78_MNI.nii.gz");
img1 = load_nii(file1{1});
aal = img1.img;

file2 = gunzip("./data/atlas/schaefer100x7MNI.nii.gz");
img2 = load_nii(file2{1});
fn = img2.img;

for i=1:9
    fn(find(fn==i))=1;
end
for i=10:15
    fn(find(fn==i))=2;
end
for i=16:23
    fn(find(fn==i))=3;
end
for i=24:30
    fn(find(fn==i))=4;
end
for i=31:33
    fn(find(fn==i))=5;
end
for i=34:37
    fn(find(fn==i))=6;
end
for i=38:50
    fn(find(fn==i))=7;
end
for i=51:58
    fn(find(fn==i))=1;
end
for i=59:66
    fn(find(fn==i))=2;
end
for i=67:73
    fn(find(fn==i))=3;
end
for i=74:78
    fn(find(fn==i))=4;
end
for i=79:80
    fn(find(fn==i))=5;
end
for i=81:89
    fn(find(fn==i))=6;
end
for i=90:100
    fn(find(fn==i))=7;
end

a=unique(aal);

N_areas=78;

% RSN projected on AAL1 cortical atlas
for n=1:78
    indn=find(aal==a(n+1));
    for Net=1:7
        FN_AAL(Net,n)=numel(find(fn(indn)==Net));
    end
end

% figure
Order=[1:2:78 78:-2:2];
for Net=1:7
    subplot(3,7,[Net Net+7])
    hold on
    barh(squeeze(FN_AAL(Net,Order)),'EdgeColor','none','Barwidth',.5,'FaceColor','r')
    ylim([0 79])
    grid on
    set(gca,'YTick',1:N_areas,'Fontsize',8)    
    set(gca,'YTickLabel',[])
end

% load cluster-centroid WEiDA eigenvectors for states
VLeida=csvread("./WEiDA4_atlasAAL78_DBS/Vemp.csv");
VLeida=VLeida.';
k=4;

Rpos=zeros(k,7);
Ppos=zeros(k,7);
Rneg=zeros(k,7);
Pneg=zeros(k,7);

for FLeida=1:k
    V=VLeida(FLeida,:);
   
    for NetYeo=1:7
        [cc p]=corr((V.*(V>0))',(FN_AAL(NetYeo,:))','type', 'Spearman');
        Rpos(FLeida,NetYeo)=cc;
        Ppos(FLeida,NetYeo)=p;
        
        [cc p]=corr((-V.*(V<0))',(FN_AAL(NetYeo,:))','type', 'Spearman');
        Rneg(FLeida,NetYeo)=cc;
        Pneg(FLeida,NetYeo)=p;
    end
end

Ppos_bf=Ppos*length(Ppos(:));
Pneg_bf=Pneg*length(Pneg(:));

csvwrite([DirOut 'OverlapRSN_Rpos.csv'],Rpos);
csvwrite([DirOut 'OverlapRSN_Rneg.csv'],Rneg);
csvwrite([DirOut 'OverlapRSN_Ppos.csv'],Ppos);
csvwrite([DirOut 'OverlapRSN_Pneg.csv'],Pneg);

% Figure of the overlap bars

figure
bar(1:k,Rpos)

figure
bar(1:k,Rneg)


