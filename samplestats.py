# -*- coding: utf-8 -*-
# samplestats.py
"""
Created on Fri Jul 13 11:57:28 2018

@author: wenyi
"""


import warnings
import math

__all__ = ['sample_stats']

#####################################
#       INFERENTIAL STATISTICS      #
#####################################

def sample_stats(x):
    """
    Compute sample statistics
    
    From a list of values to compute mean, variance with the implementation of
    a naive algorithm
    
    Parameters
    ----------
    
    Returns
    -------
    sample mean : float value
    variance : float value
    std_error : float value
    size : int
    warning : str as a warning message if not enough values
    
    References
    ----------
    
    
    Examples
    --------
    Suppose we have a list of fload values as a sample::
        362.8, 337.9, 341.4, 338.8, 285.1, 336.8, 343.0, 340.0, 339.4, 324.2
    Apply this data set to compute sample statistics:
        
    >>> from cebspy.stats.samplestats import sample_stats
    >>> sample = [362.8, 337.9, 341.4, 338.8, 285.1, 336.8, 343.0, 340.0, 339.4, 324.2]
    >>> sample_stats(sample)
    {'mean': 334.94, 'variance': 394.98488888886965, 'size': 10}
    {'mean': 334.94, 'variance': 394.98488888886965, 'standard_deviation': 19.874226749457943, 'size': 10}
    # A warning message will be returned if the sample size is less than 2
    >>> sample = []
    >>> rst = sample_stats(sample)
    >>> rst
    {'warning': 'At least two values are required for a sample'}
    >>>
    """
    (n, sum, sum_sq) = (0, 0, 0)
    if (len(x) > 1):
        for xi in x:
            n += 1
            sum += xi
            sum_sq += xi * xi
    if n > 1:
        variance = (sum_sq - (sum * sum) / n) / (n - 1)
        mean = sum / n
        return {'mean': mean, 'variance': variance, 
                'std_error':math.sqrt(variance/n), 'size': n}
    else:
        warnings.warn("At least two values are required for a sample")
    
    return {'warning':'At least two values are required for a sample'}