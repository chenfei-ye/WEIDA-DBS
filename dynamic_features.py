import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

def get_fdr(p_list):
    fdr = multipletests(p_list, method="fdr_bh", alpha=0.05)[1]
    
    return fdr

df = pd.read_csv("./data/tTIS/UPDRS3_active.csv", header=None)    # UPDRS3 socres for patients under conditions of tTIS-PRE and tTIS-ACT

# calculate symptom improvements
scl = np.zeros([11,6])                          
for i in range(11):
    scl[i,0] = -(df.iloc[i,7] - df.iloc[i,0]) / df.iloc[i,0]
    scl[i,1] = -(df.iloc[i,8] - df.iloc[i,1]) / df.iloc[i,1]
    scl[i,2] = -(df.iloc[i,9] - df.iloc[i,2]) / df.iloc[i,2]
    scl[i,3] = -(df.iloc[i,10] - df.iloc[i,3]) / df.iloc[i,3]
    scl[i,4] = -(df.iloc[i,11] - df.iloc[i,4]) / df.iloc[i,4]
    scl[i,5] = -(df.iloc[i,12] - df.iloc[i,5]) / df.iloc[i,5]
df_fea = pd.DataFrame(scl, columns=['UPDRS3','SUM','rigidity','bradykinesia','axis','tremor'])

# add state transition probability ACT-PRE changes
for i in range(4):
    for j in range(4):
        list1=[]
        list2=[]
        name= 'Ptr' + str(i+1) + 'to' + str(j+1)
        for s in range(1,12):
            pth1 = "./WEiDA4_atlasAAL78_tTIS/Ptrans" + str(s) + ".csv"        # tTIS-PRE
            data1 = pd.read_csv(pth1,header=None)
            pth2 = "./WEiDA4_atlasAAL78_tTIS/Ptrans" + str(22+s) + ".csv"     # tTIS-ACT
            data2 = pd.read_csv(pth2,header=None)    
            df_fea.loc[s-1,name]=data2.iloc[i,j] - data1.iloc[i,j]

# calculate state fractional occupancy ACT-PRE changes
data1 = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/P1emp.csv", header=None)
data2 = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/P3emp.csv", header=None)
df_fea['P1'] = data2.iloc[:,0] - data1.iloc[:,0]
df_fea['P2'] = data2.iloc[:,1] - data1.iloc[:,1]
df_fea['P3'] = data2.iloc[:,2] - data1.iloc[:,2]
df_fea['P4'] = data2.iloc[:,3] - data1.iloc[:,3]

# calculate state dwell time ACT-PRE changes
data1 = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/LT1emp.csv", header=None)
data2 = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/LT3emp.csv", header=None)
df_fea['LT1'] = data2.iloc[:,0] - data1.iloc[:,0]
df_fea['LT2'] = data2.iloc[:,1] - data1.iloc[:,1]
df_fea['LT3'] = data2.iloc[:,2] - data1.iloc[:,2]
df_fea['LT4'] = data2.iloc[:,3] - data1.iloc[:,3]

# calculate functional segregation ACT-PRE changes
data = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/network_coordination.csv", index_col=0)
df_fea['FS'] = data.loc[11:21,'functional_segregation'] - data.loc[0:10,'functional_segregation']
df_fea['wFC'] = data.loc[11:21,'within_module_fc'] - data.loc[0:10,'within_module_fc']
df_fea['aFC'] = data.loc[11:21,'across_module_fc'] - data.loc[0:10,'across_module_fc']

# calculate state metastability ACT-PRE changes
data = pd.read_csv("./WEiDA4_atlasAAL78_tTIS/Metastability_tTIS.csv", header=None)
df_fea['MS1'] = data.iloc[11:,0] - data.iloc[:11,0]    
df_fea['MS2'] = data.iloc[11:,1] - data.iloc[:11,1]    
df_fea['MS3'] = data.iloc[11:,2] - data.iloc[:11,2]    
df_fea['MS4'] = data.iloc[11:,3] - data.iloc[:11,3]   
print(df_fea) 

df_fea.to_csv("./WEiDA4_atlasAAL78_tTIS/dynamic_features.csv", index=False)

########################################### correlations #######################################################
# symptom - state dynamics
symptom = ['UPDRS3','rigidity','bradykinesia','axis','tremor']
pval = []
rho = []
for sp in symptom:
    r, p = stats.spearmanr(df_fea['Ptr1to2'],df_fea[sp])
    pval.append(p)
    rho.append(r)
pval_fdr = get_fdr(pval)

# State dynamics - MS
r_ms, p_ms = stats.spearmanr(df_fea['Ptr1to2'],df_fea['MS2'])

model=smf.ols("P2~MS2",data=df_fea)
results=model.fit()
print(results.summary())

# MS - FS
r_fs, p_fs = stats.spearmanr(df_fea['MS2'],df_fea['FS'])








    





