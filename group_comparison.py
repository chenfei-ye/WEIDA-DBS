import pandas as pd
import numpy as np
import scipy.stats
from statsmodels.stats.multitest import multipletests

def get_fdr(p_list):
    fdr = multipletests(p_list, method="fdr_bh", alpha=0.05)[1]
    
    return fdr

################################## group comparison on state fractional occupancies / dwell time ########################################## 

data0 = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/P1emp.csv", header=None)     # tTIS_PRE fractional occupancies
# data0 = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/LT1emp.csv", header=None)     # tTIS_PRE dwell time
data1 = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/P3emp.csv", header=None)     # tTIS_ACT fractional occupancies
# data1 = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/LT3emp.csv", header=None)     # tTIS_ACT dwell time


df = np.zeros([4,6])
for i in range(4):
    c0 = list(data0.iloc[:,i])
    c1 = list(data1.iloc[:,i])
    
    use = []
    mean_c0 = np.mean(c0)
    mean_c1 = np.mean(c1)
    
    std_c0 = np.std(c0)
    std_c1 = np.std(c1)
    
    t, pval = scipy.stats.ttest_rel(c1, c0)
   
    use.append(mean_c0)
    use.append(std_c0)
    use.append(mean_c1)
    use.append(std_c1)
    
    use.append(t)
    use.append(pval)
    
    df[i,:] = use
    i = i+1
df = pd.DataFrame(df)

df_marker = pd.DataFrame(df.values, index=[1,2,3,4], columns=['mean1','std1','mean2','std2','tATCvsPRE','pACTvsPRE'])
df_marker['pACTvsPRE'] = get_fdr(df_marker['pACTvsPRE'])

print(df_marker)

######################################## group comparison on state transition probabilities ########################################

for i in range(4):
    for j in range(4):
        list1=[]
        list2=[]
        name= str(i+1) + 'to' + str(j+1)
        for s in range(1,12):                   
            pth1 = "./WEiDA4_atlasAAL78_tTIS/Ptrans" + str(s) + ".csv"     # transition probability matrix for subject s under condition of tTIS-PRE
            data1 = pd.read_csv(pth1,header=None)
            list1.append(data1.iloc[i,j])
            pth2 = "./WEiDA4_atlasAAL78_tTIS/Ptrans" + str(22+s) + ".csv"   # transition probability matrix for subject s under condition of tTIS-ACT
            data2 = pd.read_csv(pth2,header=None)    
            list2.append(data2.iloc[i,j])

        tval, pval = scipy.stats.ttest_rel(list2, list1)
        print(name)
        print(np.mean(list1), np.std(list1))
        print(np.mean(list2), np.std(list2))
        print(tval, pval)

######################################### group comparison on state metastability ################################################

data = pd.read_csv("./WEiDA4_atlasAAL78_DBS/Metastability_HC_DBS.csv", header=None)

df = np.zeros([4,12])
for i in range(4):
    c0 = list(np.array(data.iloc[0:16,i]))       # HC
    c1 = list(np.array(data.iloc[16:41,i]))      # DBS-OFF
    c2 = list(np.array(data.iloc[41:66,i]))      # DBS-ON

    use = []
    mean_c0 = np.mean(c0)
    mean_c1 = np.mean(c1)
    mean_c2 = np.mean(c2)
    std_c0 = np.std(c0)
    std_c1 = np.std(c1)
    std_c2 = np.std(c2)
   
    t0, pval0 = scipy.stats.ttest_ind(c1, c0)
    t1, pval1 = scipy.stats.ttest_rel(c2, c1)
    t2, pval2 = scipy.stats.ttest_ind(c2, c0)
   
   
    use.append(mean_c0)
    use.append(std_c0)
    use.append(mean_c1)
    use.append(std_c1)
    use.append(mean_c2)
    use.append(std_c2)
    
    
    use.append(t0)
    use.append(pval0)
    use.append(t1)
    use.append(pval1)
    use.append(t2)
    use.append(pval2)
    
    
    df[i,:] = use
    i = i+1
df = pd.DataFrame(df)

df_marker = pd.DataFrame(df.values, index=[1,2,3,4], columns=['mean_HC','std_HC','mean_STNoff','std_STNoff','mean_STNon','std_STNon','tSTNoffvsHC','pSTNoffvsHC','tSTNonvsSTNoff','pSTNonvsSTNoff','tSTNonvsHC','pSTNonvsHC'])
df_marker['pSTNoffvsHC'] = get_fdr(df_marker['pSTNoffvsHC'])
df_marker['pSTNonvsSTNoff'] = get_fdr(df_marker['pSTNonvsSTNoff'])
df_marker['pSTNonvsHC'] = get_fdr(df_marker['pSTNonvsHC'])

print(df_marker)