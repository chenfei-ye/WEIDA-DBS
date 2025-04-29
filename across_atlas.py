from neuromaps.resampling import resample_images
import numpy as np
import pandas as pd
from neuromaps import stats
import nibabel as nib
from neuromaps import nulls
import copy

# WEiDA states based on AAL1 cortical atlas
state_aal = pd.read_csv("./WEiDA4_atlasAAL78_DBS/Vemp.csv", header=None)    # Cluster-centroid WEiDA eigenvectors for 4 MS states
img_org = nib.load("./data/atlas/AAL78_MNI.nii.gz")                  # AAL1 cortical atlas
img = img_org.get_fdata()
all_labels = np.unique(img)

for j in range(4):

    check_img = copy.copy(img)

    for i in range(1,79):
        check_img[img == all_labels[i]] = state_aal[j][i-1]

    imga = nib.Nifti1Image(check_img, img_org.affine)
    pth = "./WEiDA4_atlasAAL78_DBS/state" + str(j) + ".nii.gz"
    nib.save(imga, 'pth')

# WEiDA states based on schaefer100 atlas
state_schaefer = pd.read_csv("./WEiDA4_atlasSchaefer100_DBS/Vemp.csv", header=None)    # Cluster-centroid WEiDA eigenvectors for 4 MS states
img_org = nib.load("./data/atlas/schaefer100x7MNI.nii.gz")                  # AAL1 cortical atlas
img = img_org.get_fdata()
all_labels = np.unique(img)

for j in range(4):

    check_img = copy.copy(img)

    for i in range(1,101):
        check_img[img == all_labels[i]] = state_schaefer[j][i-1]

    imga = nib.Nifti1Image(check_img, img_org.affine)
    pth = "./WEiDA4_atlasSchaefer100_DBS/state" + str(j) + ".nii.gz"
    nib.save(imga, 'pth')

# topographic correlations between WEiDA states based on AAL1 and Schaefer100 atlas
state = ['state1','state2','state3','state4']
rho = np.zeros([4,1])
pvals = np.zeros([4,1])
for i in range(4):
    aal = "./WEiDA4_atlasAAL78_DBS/state" + str(j) + ".nii.gz"
    schaefer = "./WEiDA4_atlasSchaefer100_DBS/state" + str(j) + ".nii.gz"

    # transfer from MNI152 space to fsaverage10k space
    aal_fsav, schaefer_fsav = resample_images(src=aal, trg=schaefer,
                                        src_space='MNI152', trg_space='MNI152',
                                        method='linear', resampling='transform_to_alt',
                                        alt_spec=('fsaverage', '10k'))
    
    # spatial autocorrelation-preserving null model
    rotated = nulls.alexander_bloch(aal_fsav, atlas='fsaverage', density='10k',
                                n_perm=10000, seed=3512)
    df_rotated = pd.DataFrame(rotated)
    df_rotated.to_csv("./WEiDA4_atlasAAL78_DBS/Rotated3512_State" + str(i+1) + "_fsaverage10k_perm10k.csv")

    corr, pval = stats.compare_images(aal_fsav, schaefer_fsav, nulls=rotated, metric='spearmanr')
    print(f'r = {corr:.3f}, p = {pval:.3f}')
    rho[i] = corr
    pvals[i] = pval

exportData = pd.DataFrame({'state':state})
exportData["rho"] = rho 
exportData["pValue"] = pvals
exportData.to_csv("./WEiDA4_atlasAAL78_DBS/WEiDA_acrossatlas_spacorr_fsaverage10k_spearman_10kperm.csv")