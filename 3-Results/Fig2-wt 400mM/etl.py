# import library
import os, glob, sys
import numpy as np
import pandas as pd

class Data():
    DATA_DIR = '/Users/connylin/Dropbox/CA/_Publications/Manuscript RL Alcohol hab model slo1/rankinlab_slo1paper/3-Results/Fig2-WT 400mM/data' 
    MEASURES = ['RevFreq','RevSpeed','RevDur']
    pvlimit = 0.001
    alpha = 0.05

    def __init__(self):
    
    def load_data(self, **kwargs):
        csv_name = 


    def get_config():
        dir_path = os.getcwd()
        config_file = dir_path+'/config.py'
        dir_top = 'Dropbox/CA'
        while not os.path.isfile(config_file):
            # go up a level
            dir_path = os.path.dirname(dir_path)
            # break if hit dir_top
            if dir_path == dir_top:
                print(f'can not find config.py\nreaching dir_top: {dir_top}')
                break
            if dir_path == '/':
                print('can not find config.py')
                break
            # get next level config
            config_file = dir_path+'/config.py'
            print(config_file)
        if os.path.isfile(config_file):
            print(f'found config here: {config_file}')
            # import config
            sys.path.insert(0, dir_path)

    get_config()