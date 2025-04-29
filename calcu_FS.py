import pandas as pd
import numpy as np
import math
from scipy.stats import pearsonr
import glob
import scipy.stats

def calculate_functional_segregation_strength(ts, modu1, modu2):
    fc1=[]
    fc2=[]
    fc12=[]
    for i in range(len(modu1)):
        for j in range(i+1,len(modu1)):
            corr_coeff, p_value = pearsonr(ts.iloc[:,modu1[i]],ts.iloc[:,modu1[j]])
            fc1.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
    for i in range(len(modu2)):
        for j in range(i+1,len(modu2)):
            corr_coeff, p_value = pearsonr(ts.iloc[:,modu2[i]],ts.iloc[:,modu2[j]])
            fc2.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
    for i in range(len(modu1)):
        for j in range(len(modu2)):
            corr_coeff, p_value = pearsonr(ts.iloc[:,modu1[i]],ts.iloc[:,modu2[j]])
            fc12.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
    
    within_module_fc = fc1 + fc2
    across_module_fc = fc12
            
    # Calculate the average within-module FC
    average_within_module_fc = np.mean(within_module_fc)
    
    # Calculate the average across-module FC
    average_across_module_fc = np.mean(across_module_fc)
    
    # Calculate the functional segregation strength as the difference between the averages
    functional_segregation_strength = average_within_module_fc - average_across_module_fc
    
    return average_within_module_fc, average_across_module_fc, functional_segregation_strength

def calculate_nodal_functional_segregation_strength(ts, modu1, modu2, num_parcel):
    
    node_fs = [None] * num_parcel
    
    for i in range(num_parcel):
        fc_in = []
        fc_ac = []
        if i in modu1:
            for j in modu1:
                if j != i:
                    corr_coeff, p_value = pearsonr(ts.iloc[:,i],ts.iloc[:,j])
                    fc_in.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
            for j in modu2:
                corr_coeff, p_value = pearsonr(ts.iloc[:,i],ts.iloc[:,j])
                fc_ac.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
        if i in modu2:
            for j in modu1:
                corr_coeff, p_value = pearsonr(ts.iloc[:,i],ts.iloc[:,j])
                fc_ac.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
            for j in modu2:
                if j != i:
                    corr_coeff, p_value = pearsonr(ts.iloc[:,i],ts.iloc[:,j])
                    fc_in.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
                
        node_fs[i] = np.mean(fc_in) - np.mean(fc_ac)
    
    return node_fs

vemp = pd.read_csv("./WEiDA4_atlasAAL78_DBS/Vemp.csv", header=None)    # Cluster-centroid WEiDA eigenvectors for states

modu1=[]
modu2=[]
for i in range(np.shape(vemp)[0]):
    if vemp.iloc[i,1]>0:      # positive module in state 2
        modu1.append(i)
    else:
        modu2.append(i)

# BOLD time series
ti_pre = glob.glob("./data/tTIS/sub-*_active_run-01_BOLD_AAL78_MNI_basic36_sm6.csv")
ti_post = glob.glob("./data/tTIS/sub-*_active_run-03_BOLD_AAL78_MNI_basic36_sm6.csv")
ti = ti_pre + ti_post

# calculate functional segregation for each subject
functional_segregation = []
within_module_fc = []
across_module_fc=[]
for s in ti:
    ts = pd.read_csv(s,index_col=0)
    wm,am,fs = calculate_functional_segregation_strength(ts, modu1, modu2)
    within_module_fc.append(wm)
    across_module_fc.append(am)
    functional_segregation.append(fs)

tval, pval = scipy.stats.ttest_rel(functional_segregation[11:], functional_segregation[:11])     # group comparison for tTIS-ACT vs. tTIS-PRE (11 subjects)
print(np.mean(functional_segregation[:11]), np.std(functional_segregation[:11]))      # tTIS-PRE
print(np.mean(functional_segregation[11:]), np.std(functional_segregation[11:]))      # tTIS-ACT       
print(tval, pval)

net_coord = pd.DataFrame()
net_coord['within_module_fc'] = within_module_fc
net_coord['across_module_fc'] = across_module_fc
net_coord['functional_segregation'] = functional_segregation
net_coord.to_csv("./WEiDA4_atlasAAL78_tTIS/network_coordination.csv")

# calculate nodal-FS for each subject
nodal_functional_segregation = np.zeros([len(ti),np.shape(vemp)[0]])
for i in range(len(ti)):
    ts = pd.read_csv(ti[i],index_col=0)
    node_fs = calculate_nodal_functional_segregation_strength(ts, modu1, modu2, np.shape(vemp)[0])
    nodal_functional_segregation[i,:] = node_fs

nFS = pd.DataFrame(nodal_functional_segregation)
nFS.to_csv("./WEiDA4_atlasAAL78_tTIS/nodal_FS.csv")
