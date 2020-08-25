# ultrafastFitFunctions
Collection of common fit functions for discribing ultrafast and static physics phenomena to be compatible with lmfit.

To be compatible with the great fitting package lmfit, all functions need to have the independent variable as the first argument and it needs to be called 'x'.
If it is not called 'x', when calling lmfit's fit, the new name needs to be specified, which is prone to errors.


=========Usage Example ===============

    from ultrafastFitFunctions import *
    Decays.expDecay(1,1,1)

======================================




========================================================
Valid parameters in order of appearance in the functions:

    mu             :   shift of x by x0 / mu
    sigT           :   T = '', S, H, ... ; width / scaling in x-directions. Similar / related to a FWHM.
    taun           :   n = 1,2,... ; similar to sig; used as taun**1 in denominator of exp
    A,B,C,D,...    :   multiplier / normalisation    
    c,a,b,d,...    :   offsets
    k              :   slope
    other:

    alpha, beta, q : ratios

========================================================

Docstring:

    r"""
    What does the function do?

    Parameters
    ----------
    x : array_like
        `x` represents the x-coordinates of a set of datapoints.
    a : int
        `a` represents the variable, which does... .

    Returns
    -------
    y : ndarray
        An array of the same shape as x.

    """


========================================================

2Do List

1)
    Most functions are written without a mu, should be changed?
    
2)
    
    
3)
    
    
4)
