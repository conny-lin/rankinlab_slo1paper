# text report of stats results in manuscript written format
import scipy.integrate as integrate
import pandas as pd
import numpy as np

# Stats result reporting -----------------------------------------------
def print_stats_apa(stattype, dfreedom, statvalue, pvalue, **kwargs):
    # check input
    assert isinstance(stattype, str), 'stattype must be a string'
    assert isinstance(dfreedom, (int, float)), 'dfreedom must be numeric'
    if isinstance(dfreedom, float):
        dfreedom = int(dfreedom)
    assert isinstance(statvalue, (int, float)), 'statvalue must be numeric'
    # kwargs
    pvalue_limit = kwargs.pop('pvalue_limit', 0.001)
    alpha = kwargs.pop('alpha', 0.05)
    separator = kwargs.pop('separator', ' ')
    show = kwargs.pop('show', True)
    # main function
    if pvalue < pvalue_limit:
        stat_report_str = f'{stattype}({dfreedom:.0f}){separator}={separator}{statvalue:.3f}, p{separator}<{separator}{pvalue_limit}'
    elif pvalue > alpha:
        stat_report_str = f'{stattype}({dfreedom:.0f}){separator}={separator}{statvalue:.3f}, p{separator}={separator}n.s.'
    else:
        stat_report_str = f'{stattype}({dfreedom:.0f}){separator}={separator}{statvalue:.3f}, p{separator}={separator}{pvalue:.3f}'
    # print report
    if show:
        print(stat_report_str)
    return stat_report_str

def pvalue_string(pv, pvlimit=0.001, alpha=0.05, digit=3):
    if pv < pvlimit:
        pv_string = f'p<{pvlimit}'
    elif pv > alpha:
        pv_string = 'p=n.s.'
    else:
        pv = round(pv, digit)
        pv_string = f'p={pv}'
    return pv_string
# Stats result reporting END ----------------------------------------------


# Integral Class -----------------------------------------------------------
class Integral():
    def __init__(self, data):
        self.data = data

    def take_integral(self, **kwargs):
        method = kwargs.pop('method','simps')
        assert isinstance(method, str), 'method must be a string'
        if method != 'simps':
            assert False, 'this function currently does not support methods other than simps'
        result = integrate.simps(self.data)
        return result
    
    def bycolumns(self, column_names):
        output = dict()
        for c in column_names:
            output[c] = self.data[c].apply(integrate.simps, axis=1).astype('float64')
        output_df = pd.DataFrame(output)
        return output_df

    





