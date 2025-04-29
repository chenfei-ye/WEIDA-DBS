import abagen
import pandas as pd
import scipy.stats as stats
import numpy as np
import nibabel as nib
import copy
from neuromaps import nulls
from neuromaps.stats import compare_images

files = abagen.fetch_microarray(donors='all', data_dir='./Image_data/abagen-data/microarray')

# remove subcortical areas from AAL1
img_org=nib.load("./data/atlas/AAL1_MNI.nii.gz")
img = img_org.get_fdata()
all_labels = np.unique(img)

check_img = copy.copy(img)
for i in [37,38,41,42,71,72,73,74,75,76,77,78]:
    check_img[img==all_labels[i]] = 0
for i in range(91,117):
    check_img[img==all_labels[i]] = 0
    
check_all_labels = np.unique(check_img)

imga = nib.Nifti1Image(check_img, img_org.affine)
nib.save(imga, "./data/atlas/AAL78_MNI.nii.gz")

# Organize Lookuptable
org_info = pd.read_csv("./data/atlas/AAL1_MNI.csv")
idx=[36,37,40,41,70,71,72,73,74,75,76,77]
for i in range(90,116):
    idx.append(i)
info=org_info.drop(index=idx).reset_index(drop=True)

info['Index'] = info['Intensity']
info = info.drop(columns=['Space','Intensity'])
info.columns=['id','label','hemisphere']
info['structure']='cortex'

info.to_csv("./data/atlas/abagen_AAL78_Lookuptable.csv", index=False)

# gene expression projected on AAL1 cortical atlas (with interpolation by nearest centroids)
atlas_AAL78 = {'image':"./data/atlas/AAL78_MNI.nii.gz", \
                    'info':"./data/atlas/abagen_AAL78_Lookuptable.csv"}
abagen.images.check_atlas(atlas_AAL78['image'], atlas_AAL78['info'])
expression = abagen.get_expression_data(atlas_AAL78['image'], atlas_AAL78['info'], missing='centroids')

expression.to_csv("./data/atlas/AAL78_expression_centroids.csv")

# Define Transmitters into lists
acetylcholine = ["CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA6", "CHRNA7", "CHRNA10", "CHRNB1", "CHRNB2"]
dopamine = ["DRD1", "DRD2", "DRD4"]
G_aminobutyric_acid = ['GABARAP','GABARAPL1','GABARAPL2','GABARAPL3']
glutamate=['GLUD1','GLUD2','GLUL']
# Select and save transmitters only
allTransmitter = G_aminobutyric_acid + acetylcholine + dopamine + glutamate
expression = expression[allTransmitter]

# generate nodal-FS T-map aross subject groups 
nodal_fs = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/nodal_FS.csv", index_col=0)
num_parcel=78
t = [None] * num_parcel
for i in range(num_parcel):
    t[i], p = stats.ttest_rel(nodal_fs.iloc[11:22,i], nodal_fs.iloc[0:11,i])
expression["nodal_FS_T_stat"] = t

expression.to_csv("./WEiDA4_atlasAAL78_tTIS/Transmitter_nodalFS_Tmap.csv")

# topographic correlations between the nodal-FS T map and gene expression maps
fs_parc = np.array(expression["nodal_FS_T_stat"])
parcellation = nib.load('./data/atlas/AAL78_MNI.nii.gz')

# spatial autocorrelation-preserving null model
rotated = nulls.burt2020(fs_parc, atlas='MNI152', density='2mm',
                                n_perm=10000, seed=3512, parcellation=parcellation)
df_rotated = pd.DataFrame(rotated)
df_rotated.to_csv("./WEiDA4_atlasAAL78_tTIS/Rotated3512_MNIparc78_perm10k_FStmap.csv")

rho = np.zeros(len(allTransmitter))
pvals = np.zeros(len(allTransmitter))
for j in range(len(allTransmitter)):
    gene_parc = np.array(expression[allTransmitter[j]])
    corr, pval = compare_images(fs_parc, gene_parc, nulls=rotated, metric='spearmanr')
    print(f'r = {corr:.3f}, p = {pval:.3f}')
    rho[j] = corr
    pvals[j] = pval

exportData = pd.DataFrame({'Gene':allTransmitter})
exportData["rho"] = rho 
exportData["pValue"] = pvals
exportData.to_csv("./WEiDA4_atlasAAL78_tTIS/Transmitter_nodalFS_spacorr_MNIparc78_perm10k_spearmanr.csv")


