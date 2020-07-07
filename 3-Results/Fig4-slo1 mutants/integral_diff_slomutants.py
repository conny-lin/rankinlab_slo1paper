# integral anova
# Conny | July 7, 2020

# import library
print('importing libraries')
import os, glob, sys
import numpy as np
import pandas as pd
import scipy.integrate as integrate
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats import multicomp
import matplotlib.pyplot as plt
import seaborn as sns
# add local library
lib = '/Users/connylin/Dropbox/CA/_Publications/Manuscript RL Alcohol hab model slo1/rankinlab_slo1paper/2-Materials-Methods/local_lib'
if lib not in sys.path:
    sys.path.insert(1, lib)
import etl, stats, graphpack

# set variables
print('set variables')
DATA_DIR = '/Users/connylin/Dropbox/CA/_Publications/Manuscript RL Alcohol hab model slo1/rankinlab_slo1paper/3-Results/Fig4-slo1 mutants/data'
OUTPUT_DIR_NAME = 'Score_Integral'
strain_dirs = glob.glob(DATA_DIR+'/*')
sub_dir = 'TWR/Dance_ShaneSpark4'
PRJ_TAG = 'integral_diff' 
MEASURES = ['RevFreq','RevSpeed','RevDur']
pvlimit = 0.001
alpha = 0.05
independent_variable = 'groupname'

for strain_dir in strain_dirs:
    strain_name = os.path.basename(strain_dir)
    print(f'\n\nProcessing: {strain_name}')
    # data directory
    raw_data_dir = os.path.join(strain_dir, sub_dir)
    # make output folder
    output_dir = os.path.join(os.path.dirname(raw_data_dir),
                                OUTPUT_DIR_NAME)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    # get data
    rawdata, db = etl.get_data(raw_data_dir)
    df_transform = rawdata.pivot(index='mwtid',
                                columns='tap',
                                values=MEASURES)
    # calculate integral (requires stats and etl package)
    intobj = stats.Integral(df_transform)
    data_integral = intobj.bycolumns(MEASURES)
    data_integral = etl.merge_data_mwtdb(data_integral, db)

    # save excel graphing output by measures---
    iv = 'groupname'
    savefname = f'graph_data_{PRJ_TAG}.xlsx'
    savepath = os.path.join(output_dir, savefname)
    graphpack.save_excel_graphdata(data_integral, 'groupname', MEASURES, savepath)

    # set up anova report ---
    filename = os.path.join(output_dir, 'anova.txt')
    text_file = open(filename, 'w')
    text_file.write('ANOVA Integral\n')
    for msr in MEASURES:
        # write header
        write_string = f'<{msr}>'
        print(write_string)
        text_file.write(write_string)
        write_string = '--ANOVA--'
        print(write_string)
        text_file.write(write_string)
        # run anova 
        ols_equation = f'{msr} ~ strain + rx + strain:rx'
        lm = ols(ols_equation, data_integral).fit()
        anovaT = anova_lm(lm)
        # print ANOVA results
        for factor in anovaT.index[:-1]:
            df = anovaT.loc[factor, 'df']
            F = anovaT.loc[factor, 'F']
            pv = stats.pvalue_string(anovaT.loc[factor, 'PR(>F)'], 
                                    pvlimit, alpha)
            write_string = f'{factor}: F({df:.0f})={F:.3f}, {pv}'
            print(f'\t{write_string}')
            text_file.write(write_string)
        # posthoc (get rid of nan)
        data_no_nan = data_integral[[independent_variable, msr]].dropna()
        r = multicomp.pairwise_tukeyhsd(data_no_nan[msr], 
                                        data_no_nan[independent_variable])
        # print posthoc results
        write_string = '--Tukey HSD posthoc--'
        print(write_string)
        text_file.write(write_string)
        df = pd.DataFrame(data=r._results_table.data[1:], 
                            columns=r._results_table.data[0])
        for i in df.index.values:
            g1 = df.loc[i, 'group1']
            g2 = df.loc[i, 'group2']
            pv = stats.pvalue_string(df.loc[i, 'p-adj'], pvlimit, alpha)
            write_string = f'\t{g1} vs {g2}, {pv}'
            print(f'\t{write_string}')
            text_file.write(write_string)
        
        # swarm plot --
        # make graph
        plt.figure()
        sns.swarmplot(x='groupname',y=msr,data=data_integral)
        plt.title(msr)
        # construct save path
        savefname = f'{msr}_swarmplot_{PRJ_TAG}.jpeg'
        savepath = os.path.join(output_dir, savefname)
        # save fig
        plt.savefig(savepath)
        plt.show() 
        plt.close()
    # finish text export
    text_file.close()
    print(f'\n{strain_name} complete')
    print('=============================================')

