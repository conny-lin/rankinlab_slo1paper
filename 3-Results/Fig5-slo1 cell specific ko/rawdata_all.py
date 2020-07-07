# import library
import os, glob, sys
import numpy as np
import pandas as pd

DATA_DIR = '/Users/connylin/Dropbox/CA/_Publications/Manuscript RL Alcohol hab model slo1/rankinlab_slo1paper/3-Results/Fig5-slo1 cell specific ko/data'

# glob for strain folders
strain_dirs = glob.glob(DATA_DIR+'/*')
sub_dir = 'Dance_ShaneSpark4'

# glob for strain folders
strain_dirs = glob.glob(DATA_DIR+'/*/')
sub_dir = 'Dance_ShaneSpark4'
db = []
for strain_dir in strain_dirs:
    data_path = os.path.join(strain_dir, sub_dir, 'mwtdb.csv')
    df = pd.read_csv(data_path)
    df['strain_dir'] = np.tile(strain_dir, df.shape[0])
    db.append(df)
db = pd.concat(db, axis=0, ignore_index=True)
db.to_csv(os.path.join(DATA_DIR, 'mwtdb_all.csv'), index=False)
print('database stored')


# combine raw data
data = []
for strain_dir in strain_dirs:
    # get strain data dir
    raw_data_dir = os.path.join(strain_dir, sub_dir)
    # get raw data
    rawdata = pd.read_csv(os.path.join(raw_data_dir, 'rawdata.csv'))
    # get database
    db = pd.read_csv(os.path.join(raw_data_dir, 'mwtdb.csv')) 
    # fill in NaN rx with 0mM
    db.loc[ db['rx'].isna(), 'rx'] = ['0mM']
    # merge database with raw data
    df = rawdata.merge(db[['mwtid','mwtpath','strain','rx','groupname']], 
                                left_on='mwtid', 
                                right_on='mwtid', 
                                how='left')
    df.drop(columns='mwtid', inplace=True)
    data.append(df)
# concat
data = pd.concat(data, axis=0, ignore_index=True)
# drop duplicates
data.drop_duplicates(inplace=True)
# save combined data
data.to_csv(os.path.join(DATA_DIR, 'rawdata_all.csv'), index=False)