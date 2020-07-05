# import library
import os, glob, sys
import numpy as np
import pandas as pd

# global variables
# project directory data is located in the same directory of this code
DATA_DIR = os.path.join(os.path.dirname( os.getcwd() ), 'data') 
# define measures, statistical setting
MEASURES = ['RevFreq','RevSpeed','RevDur']
pvlimit = 0.001
alpha = 0.05