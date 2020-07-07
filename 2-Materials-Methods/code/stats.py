# text report of stats results in manuscript written format
import scipy.integrate as integrate

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


def take_integral(d, **kwargs):
    method = kwargs.pop('method','simps')
    assert isinstance(method, str), 'method must be a string'
    if method != 'simps':
        assert False, 'this function currently does not support methods other than simps'
    result = integrate.simps(d)
    return result



