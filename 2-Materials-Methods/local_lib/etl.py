import os
import pandas as pd
import numpy as np

def confirm_import():
    print('confirm')

    
def get_data(raw_data_dir):
    # get raw data ---
    # get strain data dir
    # get raw data
    rawdata = pd.read_csv(os.path.join(raw_data_dir, 'rawdata.csv'))
    # get database
    db = pd.read_csv(os.path.join(raw_data_dir, 'mwtdb.csv')) 
    # fill db rx empty to 0mM
    db['rx'].fillna('0mM', inplace=True)
    return rawdata, db

def merge_data_mwtdb(data, mwtdb, **kwargs):
    columns = kwargs.pop('columns', ['mwtid','groupname', 'strain','rx'])
    on = kwargs.pop('on','mwtid')
    merged_data = mwtdb[columns].merge(data, left_on=on, right_on=on, how='right')
    return merged_data

