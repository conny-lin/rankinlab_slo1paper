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