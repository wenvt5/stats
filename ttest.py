# -*- coding: utf-8 -*-
# ttest.py
"""
A collection of basic statistial functions for CEBS.
Created on Mon Jul 2 15:38:31 2018

@author: wenyi

"""
import warnings
import math
import numpy as np
from scipy.stats import distributions

from cebspy.stats.commons import valid_floats as vfloats
from cebspy.stats.samplestats import sample_stats


__all__ = ['one_sample_t_test','t_test']

#####################################
#       INFERENTIAL STATISTICS      #
#####################################

    
def one_sample_t_test(x, popmean=0):
    """
    Compute the t statistics for one sample
    
    Calculate t-statistics of one sample comparing with a population mean
    
    Parameters
    ----------
    x : a list of numeric values as a sample
    popmean : a float/numeric value as the population mean
    

    Returns
    -------
    name : a str value as 'One Sample t-test'
    t : t statistic value
    df : degree of freedom
    p_value : probability (p) value
    
    Notes
    -----
    
    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Student%27s_t-test

    Examples
    --------
    >>> from cebspy.stats.ttest import one_sample_t_test
    >>> sample = [362.8, 337.9, 341.4, 338.8, 285.1, 336.8, 343.0, 340.0, 339.4]
    >>> results = one_sample_t_test(sample, popmean=335)
    >>> print(json.dumps(results, indent=4))
    {
        "name": "One Sample t-test",
        "sample_stats": {
            "mean": 336.1333333333334,
            "variance": 428.3374999999796,
            "std_error": 6.89877203243833,
            "size": 9
        },
        "t": 0.1642804441144597,
        "df": 8,
        "p_value": 0.8735851934842374,
        "std_error": 6.89877203243833
    }
    >>>

    """
    x_stats = sample_stats(x)
    n = x_stats.get('size')
    results = {'name':'One Sample t-test',
               'sample_stats':x_stats,}
    df = n - 1
    if df > 0:
        mean = x_stats.get('mean')
        std_error = x_stats.get('std_error')
        if std_error > 0:
            t = (mean - popmean) / std_error
            prob = distributions.t.sf(np.abs(t), df) * 2
            results.update({'t':t,
                            'df':df,
                            'p_value':prob,
                            'popmean':popmean,
                            'std_error':x_stats.get('std_error')})
        else:
            results['warnings'] = 'No variance: data are essentially constant'
    else:
        results['warnings'] = 'The sample size need to be at least 2'
    return results
    

def equal_var_t_test_estimates(variance1, n1, variance2, n2):
    """
    Compute the standard error based on pooled variance
    
    Calucate pooled variance and standard error based on the two samples sizes
    assumming variances are equal
    return the standard errors and the degree of freedom
    """
    df = n1 + n2 - 2
    sum_variance = (n1 - 1) * variance1 + (n2 - 1) * variance2
    sp = math.sqrt(sum_variance / df)
    std_error = sp * math.sqrt(1 / n1 + 1 / n2)
    return (df, std_error);

def unequal_var_t_test_estimates(variance1, n1, variance2, n2):
    """
    Compute the standard error based on pooled variance
    
    Calucate pooled variance and standard error based on the two samples sizes
    assuming variances are not equal
    return the pooled standard errors and the pooled degree of freedom
    """
    denom = (variance1 / n1) ** 2 / (n1 - 1) + (variance2 / n2) ** 2 / (n2 - 1)
    pooled_variance = variance1 / n1 + variance2 / n2
    df = pooled_variance ** 2 / denom
    std_error = math.sqrt(pooled_variance)
    return(df, std_error)


def two_samples_t_test(x, y, var_equal):
    """
    Compute t-test statistics for two samples
    """
    x_stats = sample_stats(x)
    y_stats = sample_stats(y)
    name = 'Two Sample t-test'
    if var_equal:
        (df, std_error) = equal_var_t_test_estimates(
                x_stats.get('variance'), x_stats.get('size'),
                y_stats.get('variance'), y_stats.get('size'))
    else:
        (df, std_error) = unequal_var_t_test_estimates(
                x_stats.get('variance'), x_stats.get('size'),
                y_stats.get('variance'), y_stats.get('size'))
        name = 'Welch Two Sample t-test'
        
    t = (x_stats.get('mean') - y_stats.get('mean')) / std_error
    prob = distributions.t.sf(np.abs(t), df) * 2
    return {'name':name,
            'sample1_stats':x_stats, 'sample2_stats':y_stats,
            't':t, 'df':df, 'p_value':prob, 
            'std_error':std_error, 'var_equal':var_equal}


def t_test(x, y = None, alternative = 'two.sided', popmean = 0, paired = False,
           var_equal = False, conf_level = 0.95, nan_policy = 'omit'):
    """
    Performs t-test of one sample or two samples with data.
    
    This is an implementation of t-test with different cases including
    two-sided or equal, left-side or less, right-side or greater with
    variance equal as True or False in case of two samples.
    
    Parameters
    ----------
    
    x : a list of numeric data values
    
    y : a list of numeric data values or an empty list or None, optional
        
    alternative : a character string specifying the alternative hypothesis, optional
                  If defined, must be one of "two.sided" (default), "greater" or "less".
                  You can specify just the initial letter.
        Depedent on the null hypothesis of the test as
        1. Two means from the two idenpendent samples are equal
        2. The mean from the sample 1 is less than the mean from the sample two
        3. The mean from the sample 1 is greater than the mean from the sample two
    
    popmean : a number indicating the true value of the mean 
         (or difference in means if you are performing a two sample test).
         
    paired : a logical indicating whether you want a paired t-test.

    var.equal : a logical variable indicating whether to treat the two variances
                as being equal. If True then the pooled variance is used to 
                estimate the variance otherwise (False) the Welch (or Satterthwaite) 
                approximation to the degrees of freedom is used.

    conf.level : confidence level of the interval.

    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.
        
        
    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Student%27s_t-test

    Examples
    --------
    Suppose we have the summary data for two samples, as follows::
        sample1 : 362.8, 337.9, 341.4, 338.8, 285.1, 336.8, 343.0, 340.0
        sample2 : 422.2, 454.3, 429.8, 376.8, 408.7, 488.9, 405.9, 444.2, 369.1, 441.7
    Apply the two samples to t-test assume variance either equal or not (default)

    >>> from cebspy.stats.ttest import t_test
    >>> sample1 = [362.8, 337.9, 341.4, 338.8, 285.1, 336.8, 343.0, 340.0]
    >>> sample2 = [422.2, 454.3, 429.8, 376.8, 408.7, 488.9, 405.9, 444.2, 369.1, 441.7]
    >>> results = t_test(sample1, sample2) # assume var_equal = False (default)
    >>> results.keys()
    dict_keys(['method', 'has_output', 'has_error', 'output'])
    >>> import json
    >>> if 'output' in results:
    ...     print(json.dumps(results.get('output'), indent=4))
    ...
    {
        "name": "Welch Two Sample t-test",
        "sample1_stats": {
            "mean": 335.725,
            "variance": 487.8135714285475,
            "std_error": 7.808757675108662,
            "size": 8
        },
        "sample2_stats": {
            "mean": 424.15999999999997,
            "variance": 1299.0671111111685,
            "std_error": 11.397662528392251,
            "size": 10
        },
        "t": -6.400885972084148,
        "df": 15.142437005481307,
        "p_value": 1.140116336879273e-05,
        "std_error": 13.816056149990319,
        "var_equal": false
    }

    >>> results2 = t_test(sample1, sample2, var_equal=True) # variance equal
    >>> if 'output' in results2:
    ...     print(json.dumps(results2, indent=4))
    ...
    {
        "method": "T-test",
        "has_output": true,
        "has_error": false,
        "output": {
            "name": "Two Sample t-test",
            "sample1_stats": {
                "mean": 335.725,
                "variance": 487.8135714285475,
                "std_error": 7.808757675108662,
                "size": 8
            },
            "sample2_stats": {
                "mean": 424.15999999999997,
                "variance": 1299.0671111111685,
                "std_error": 11.397662528392251,
                "size": 10
            },
            "t": -6.067557130079347,
            "df": 16,
            "p_value": 1.6306979762711614e-05,
            "std_error": 14.575058479728472,
            "var_equal": true
        }
    }
    >>>
    """
    s1, cn1, ce1 = vfloats(x)
    results = {'method':'T-test', 'has_output':bool(0), 'has_errors':bool(0)}
    test = None
    
    if (len(s1) > 1):
        if y is None:
            # One sample t-test
            test = one_sample_t_test(s1, popmean)
        else:
            s2, cn2, ce2 = vfloats(y)
            if (len(s2) > 1):
                test = two_samples_t_test(s1, s2, var_equal)
            else:
                test = one_sample_t_test(s1, popmean)
            results['is_finished'] = bool(1)
    else:
        message = 'At least one sample with a minimum of two values is required'
        warnings.warn(message)
        results['has_error'] = bool(1)
        results['warnings'] = [message]

    if (test is not None):
        results['has_output'] = bool(1)
        results['output'] = test

    return results



