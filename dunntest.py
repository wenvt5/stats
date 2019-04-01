# -*- coding: utf-8 -*-
# dunntest.py
"""
Created on Thu Jul 26 15:58:00 2018

@author: wenyi
"""

import numpy as np
import pandas as pd
import warnings

import cebspy.stats._ties_utils as tu
from scipy.stats import norm

__all__ = ['dunn_test']


def dunn_test(doses, responses):
    """
    Dunn's multiple comsprison test
    
    Deternine the lowest dose different from control
    
    Parameters
    ----------
    doses : a list of folat values as doses ordered ascendingly
    
    responses : a list of float values as responses
    
    Examples
    --------
    >>> from cebspy.stats.dunntest import dunn_test as dunn
    >>> import cebspy.tests.examples as examples
    >>> doses, responses = examples.dunn_data()
    >>> if ('doses' in data and 'responses' in data):
    ...     doses = data['doses']
    ...     responses = data['responses']
    >>> if (len(doses) == len(responses)):
    ...     results = dunn(doses, responses)
    ...     if 'output' in results:
    ...         print(results)
    ...
    {'dose': [0, 1000, 5000, 10000],
     'count': [40, 39, 37, 40], 
     'dose_rank': [0, 1, 2, 3], 
     'rank_mean': [77.7, 89.41025641025641, 76.78378378378379, 70.25],
     'dunnsign': [0, 0, -1, -1],
     'mult_comp_signif': [0, 0, 0, 0]}
    
    """
    results = {'method':"Dunn's test",
               'has_output':bool(0),
               'has_errors':bool(0)}
    warn_message = None
    if (len(doses) == len(responses)):
        dose_groups = np.unique(doses) # also sorted
        if (0 in dose_groups):
            if (len(dose_groups) > 1):
                ranks = np.unique(responses, return_counts=True)
                # np.shape(ranks)
                ties = np.unique(ranks[1], return_counts=True)
                # tie adjustment correction
                correction = tu._ties_correction(ties)
                weights = tu._ranks_weights(ranks)
                #ranks2 = tuple([ranks[0], ranks[1], np.array(weights)]) # ['value', 'count', 'rank']
                ranks2 = pd.DataFrame({'values':ranks[0].tolist(),
                                       'counts':ranks[1].tolist(),
                                       'ranks':weights})
                ranks2.set_index('values', inplace=True)
                ex_df = pd.DataFrame({'doses':doses, 'response':responses})
                ex_df['counts'] = [ranks2['counts'][value] for value in responses]
                ex_df['ranks'] = [ranks2['ranks'][value] for value in responses]
                #rank_sums = ex_df.groupby('doses')['ranks'].agg(['sum','count'])
                rank_means = ex_df.groupby('doses')['ranks'].agg(['mean','count'])
                rank_means['dose_count'] = list(range(len(rank_means)))
                rank_means.reset_index(drop=False, inplace=True)
                # find variance
                n_total = len(ex_df)
                v = (n_total * (n_total + 1))/12
                # get crit values ... warnings are for CONTROL group, which is fine
                prob05 = 1 - (.05 / (2 * max(rank_means['dose_count'])))
                prob01 = 1 - (.01 / (2 * max(rank_means['dose_count'])))
                z_score05 = norm.ppf(prob05)
                z_score01 = norm.ppf(prob01)
                # loop through doses, determine significance
                dunnsigns = [0]
                mult_comp_signifs = [0]
                for i in range(1, len(rank_means), 1):
                    sig = - 1
                    if (rank_means['mean'][i] >= rank_means['mean'][0]):
                        sig = 0
                    dunnsigns.append(sig)
                    #num = rank_means['count'][i]
                    rankdiff = abs(rank_means['mean'][i] - rank_means['mean'][0])
                    comp2 = v * (1/rank_means['count'][i] + 1/rank_means['count'][0])
                    comp2 = (comp2 * (1 - correction / (n_total ** 3 - n_total))) ** .5
                    sig05 = rankdiff - (z_score05 * comp2)
                    sig01 = rankdiff - (z_score01 * comp2)
                    signific = 0
                    if (sig05 > 0):
                        signific = 1
                    if (sig01 > 0):
                        signific = 2
                    mult_comp_signifs.append(signific)
                tests = {'is_finished':bool(1),
                        'has_output':bool(1),
                        'output':{'dose':rank_means['doses'].tolist(),
                                  'count':rank_means['count'].tolist(),
                                  'dose_rank':rank_means['dose_count'].tolist(),
                                  'rank_mean':rank_means['mean'].tolist(),
                                  'dunnsign':dunnsigns,
                                  'mult_comp_signif':mult_comp_signifs}}
                results.update(tests)
            else:
                warn_message = 'No enough treatment groups'
                warnings.warn(warn_message)
        else:
            warn_message = 'The control group (dose = 0) is missing'
            warnings.warn(warn_message)
    else:
        warn_message = 'The number of values betwee doses and responses are not equal'
        warnings.warn(warn_message)
    if (warn_message is not None):
        results.update({'has_errors':bool(1),
                        'warnings':[warn_message]})

    return results
