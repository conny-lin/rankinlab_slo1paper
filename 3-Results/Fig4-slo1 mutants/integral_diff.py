import os, glob, sys
import numpy as np
import pandas as pd
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# global variables
# project directory data is located in the same directory of this code
DATA_DIR = os.getcwd() # used to be proj_dir
# define save suffix for this project
PRJ_TAG = 'integral_diff' # used to be proj_suffix
# define measures, statistical setting
MEASURES = ['RevFreq','RevSpeed','RevDur']
pvlimit = 0.001
alpha = 0.05

p = 

# get raw data
def get_data(raw_data_dir):
    # get raw data
    rawdata = pd.read_csv(os.path.join(raw_data_dir, 'rawdata.csv'))
    # get database
    db = pd.read_csv(os.path.join(raw_data_dir, 'mwtdb.csv')) 
    # fill in NaN rx with 0mM
    db.loc[ db['rx'].isna(), 'rx'] = ['0mM']

    # merge database with raw data
    df_merge = rawdata.merge(db[['mwtid','strain','rx']], 
                            left_on='mwtid', 
                            right_on='mwtid', 
                            how='left'
                            )
    return df_merge

# define groups (0mM vs 400mM) 
rxs = df_merge['rx'].unique().tolist()

# create output array
df_output = pd.DataFrame(df_db[['mwtid','rx','strain','groupname']])
for msr in MEASURES:
    df_output[msr] = 'NaN'
#set index
df_output.set_index('mwtid', inplace=True)

# transform each plate's data in rows, columns = taps
for msr in MEASURES: # for each measure
    # transform data
    df_pivot = df_merge.pivot(index='mwtid',columns='tap',values=msr)
    for i in df_pivot.index: # go through each plate
        d = df_pivot.loc[i,:]
        df_output.loc[i,msr] = integrate.simps(d)
# remove missing values
# there are some nan values. probably don't have all the data, check which ones
print(f'database has {db.shape[0]} plates')
print(f'result has {df_pivot.shape[0]} plates')
print('remove missing plates')

# remove missing values (not sure why dropna() doesn't work)
bol_output = (df_output == 'NaN')
ind_to_drop = df_output.index[bol_output.any(axis=1)]
df_output.drop(index=ind_to_drop,inplace=True)
print(f'len of output {len(df_output)}') # not sure why using dropna() didn't drop

# ensure data are in the right datatype
for msr in MEASURES: # for each measure
    df_output[msr] = df_output[msr].astype('float64').copy()

# create excel graphing output
iv = 'rx'
savefname = f'Desc_Stats_{PRJ_TAG}.xlsx'
savepath = os.path.join(DATA_DIR,savefname)
# save to excel
with pd.ExcelWriter(savepath) as writer:
    for msr in MEASURES:
        # get data
        df = df_output[[iv,msr]].copy()
        # create summary
        t = df.groupby(iv).agg(['count','mean','sem'])
        # save to excel
        t.to_excel(writer, sheet_name=msr)
    print('done')

# develop written report 
from stats.report import print_stats_apa
filename = os.path.join(DATA_DIR, 'anova_txt_report.txt')
text_file = open(filename, 'w')
text_file.write('Integral differences between 0mM and 400mM\n')
for msr in MEASURES:
    print(f'\n{msr} ----------------')
    # get data
    df = df_output[['rx',msr]].copy()
    # run anova 
    lm = ols(f'{msr} ~ rx', data=df).fit()
    anovaT = anova_lm(lm)
    display(anovaT)
    # report result in APA format
    anova_apa = print_stats_apa(anovaT.columns[3], 
                                 anovaT.loc[iv,'df'],
                                 anovaT.loc[iv,anovaT.columns[3]], 
                                 anovaT.loc[iv,'PR(>F)'], 
                                 pvalue_limit = pvlimit, 
                                 alpha=alpha,
                                 separator='',
                                 show=True)
    output_str = f'{msr}: {anova_apa}\n'
    print(output_str)
    text_file.write(output_str)
text_file.close()