import scipy.integrate as integrate

def take_integral(d, **kwargs):
    method = kwargs.pop('method','simps')
    assert isinstance(method, str), 'method must be a string'
    if method != 'simps':
        assert False, 'this function currently does not support methods other than simps'
    result = integrate.simps(d)
    return result

