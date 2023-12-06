# -*- coding: utf-8 -*-
"""
import random 
import numpy as np
import pandas as pd

def roll(a,b):
    inter_list = []
    for i in range(1000):
        a_s = np.array([])
        b_s = np.array([])
        for i in range(a):
            a_s = np.append(a_s,np.random.randint(0,29799))
        for i in range(b):
            b_s = np.append(b_s,np.random.randint(0,29799))
        inter=len(np.intersect1d(a_s,b_s))
        inter_list.append(inter)
    
    inter_df = pd.DataFrame(inter_list)
    inter_df = inter_df.sort_values(by=[0])
    return inter_df


test = pd.read_csv(r'<SET_PATH>\fst_DEG_stat_tissue.txt',sep='\t')
test = pd.read_csv(r'<SET_PATH>\fst_DCG_stat_tissue.txt',sep='\t')

test['p<0.001'] = 0
test['p<0.01'] = 0
test['p<0.05'] = 0
inter_df_all = pd.DataFrame(np.array([0]*len(test)*1000).reshape(1000,len(test)))
for i in range(len(test)):
    Fst = test.iloc[i,1]
    DEG = test.iloc[i,2]
    overlap = test.iloc[i,3]
    inter_df = roll(Fst, DEG)
    inter_df_0001 = inter_df.max()[0]
    inter_df_001 = inter_df[-10:].iloc[0][0]
    inter_df_005 = inter_df[-50:].iloc[0][0]
    if overlap > inter_df_0001:
        test.iloc[i,4] = 1
    elif overlap > inter_df_001:
        test.iloc[i,5] = 1
    elif overlap > inter_df_005:
        test.iloc[i,6] = 1
    inter_df = inter_df.reset_index(drop=True)
    inter_df_all.iloc[:,i] = inter_df
    
inter_df_all.to_csv('<SET_PATH>\E_test.txt',sep='\t',index=False)
