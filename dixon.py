# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 12:23:01 2018

@author: yingw
""" 
#os.chdir('C:\\Users\\yingw\\AnacondaProjects\\cebstats')
#import pandas as pd
import numpy as np
import qdixon

#x = [0.2022, 0.2111, 0.2173, 0.2190, 0.2268, 0.2270, 0.2334, 0.2338, 0.2338, 0.2354, 0.2371, 0.2372, 0.2378, 0.2418, 0.2451, 0.2455, 0.2460, 0.2549, 0.2550, 0.2633, 0.2644, 0.2724, 0.2915]

__all__ = ['dixon']

def dixon(x, type = 0, opposite = False, two_sided = True):
    x.sort()
    n = len(x)
    if ((type == 10 or type == 0) and (n < 3 or n > 30)):
        return("Sample size must be in range 3-30")
    if ((type == 11) and (n < 4 or n > 30)):
        return("Sample size must be in range 4-30")
    if ((type == 12) and (n < 5 or n > 30)):
        return("Sample size must be in range 5-30")
    if ((type == 20) and (n < 4 or n > 30)):
        return("Sample size must be in range 4-30")
    if ((type == 21) and (n < 5 or n > 30)):
        return("Sample size must be in range 5-30")
    if ((type == 22) and (n < 6 or n > 30)):
        return("Sample size must be in range 6-30")
    if type not in  [0, 10, 11, 12, 20, 21, 22]:
        return("Incorrect type")
    if type == 0:
        if n <= 7:
            type = 10
        elif (n > 7 and n <= 10):
            type = 11
        elif (n > 10 and n <= 13):
            type = 21
        else:
            type = 22
    if ((x[n-1] - np.mean(x)) < (np.mean(x) - x[0]) != opposite ):
        alt = 'lowest value ' + str(x[0]) +' is an outlier'
        if type == 10:
            Q = (x[1] - x[0]) / (x[n-1] - x[0])
        elif type == 11:
            Q = (x[1] - x[0]) / (x[n-2] - x[0])
        elif type == 12:
            Q = (x[1] - x[0]) / (x[n-3] - x[0])
        elif type == 20:
            Q = (x[2] - x[0]) / (x[n-1] - x[0])
        elif type == 21:
            Q = (x[2] - x[0]) / (x[n-2] - x[0])
        else:
            Q = (x[2] - x[0]) / (x[n-3] - x[0])
    else:
        alt = 'highest value ' + str(x[n-1]) + ' is an outlier'
        if type == 10:
            Q = (x[n-1] - x[n-2]) / (x[n-1] - x[0])
        elif type == 11:
            Q = (x[n-1] - x[n-2]) / (x[n-1] - x[1])
        elif type == 12:
            Q = (x[n-1] - x[n-2]) / (x[n-1] - x[2])
        elif type == 20:
            Q = (x[n-1] - x[n-3]) / (x[n-1] - x[0])
        elif type == 21:
            Q = (x[n-1] - x[n-3]) / (x[n-1] - x[1])
        else:
            Q = (x[n-1] - x[n-3]) / (x[n-1] - x[2])
    pval = qdixon.qdixon(Q, n, type)
    if two_sided:
        pval = 2 * pval
        if pval > 1:
            pval = 2 - pval
    Q_list = ['Q',Q]
    RVAL = {'statistic' : Q_list, 'alternative' : alt, 'p.value' : pval,
        'method' :"Dixon test for outliers",  'class': 'htest'}
    return(RVAL)
