# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 14:52:22 2018

@author: wenyi
"""
import numpy as np
import pandas as pd
import warnings

__all__ = ['shirley_test']

def _mult_comparison(dose_count, test_nums, test_stats, test_doses):
    ## add SAS crit values
    C01 = [0, 2.575, 2.607, 2.615, 2.618, 2.620, 2.621, 2.622]
    C05 = [0, 1.96, 2.015, 2.032, 2.040, 2.044, 2.047, 2.0485]
    
    B01 = [0, 0, 3, 4, 4, 4, 4, 4]
    B05 = [0, 0, 3, 4, 5, 6, 6, 6]
    mult_comp_signif = [None] * len(dose_count)
    nonsignif_flag = 'NO'
    if (len(dose_count) <= 7):
        for i in range(len(dose_count)):
            if (nonsignif_flag == 'NO'):
                dosenum = dose_count[i] + 1
                crit05 = C05[dosenum] - ( (B05[dosenum] / 100) * (1 - (test_nums[i] / test_nums[-1])))
                crit01 = C01[dosenum] - ( (B01[dosenum] / 100) * (1 - (test_nums[i] / test_nums[-1])))	
                if(test_stats[i] >= crit01):
                    mult_comp_signif[i] = 2
                else:
                    if (test_stats[i] >= crit05):
                        mult_comp_signif[i] = 1
                    else:
                        mult_comp_signif[i] = 0
                        nonsignif_flag = 'YES'
        return {'is_finished':bool(1),
                'has_output':bool(1),
                'output':{'dose':test_doses,
                          'mult_comp_signif':mult_comp_signif}}
    else:
        warn_message = ('Number of dose groups exceeds the miximum '
                        'number of critical values from SAS')
        warnings.warn(warn_message)
        return {'has_errors':bool(1),
                'warnings':[warn_message]}


def shirley_test(doses, responses, tau):
    """
    Shirley's doses and responses test
    
    Deternine the lowest dose different from control
    
    Parameters
    ----------
    doses : a list of folat values as doses ordered ascendingly
    
    responses : a list of float values as responses
    
    tau : a float value as the tau statistic from Kendall's correlation test
    
    References
    ----------
    
    ..[1] Shirley.R as an R script from CEBS table project
    
    Examples
    --------
    
    
    """
    if (len(doses) == len(responses)):
        results = {'method':"Shirley's test",
                   'has_output':bool(0),
                   'has_errors':bool(0)}
        dose_groups = np.unique(doses) # also sorted
        if (len(dose_groups) > 1 and 0 in dose_groups):
            #dose_ties = np.unique(doses, return_counts=True)
            #last_index = len(responses)
            groups = len(dose_groups)
            test_stats = []
            dose_count = []
            test_doses = []
            test_nums = []
            for g in range(groups - 2, -1, -1):
                max_dose = dose_groups[g + 1]
                tmp_doses = [doses[i] for i in range(len(doses)) if doses[i] <= max_dose]
                tmp_values = [responses[i] for i in range(len(doses)) if doses[i] <= max_dose]
                doses_ranks = np.unique(tmp_doses, return_counts=True)
                ranks = np.unique(tmp_values, return_counts=True)
                ties = np.unique(ranks[1], return_counts=True)
                n_total = len(tmp_values)
                correction = 0
                if (len(ties) > 1 and len(ties[1]) > 1):
                    #remove_singleton as ties number is 1
                    indices = [i for i in range(len(ties[0])) if ties[0][i] != 1]
                    ties_numbers = [ties[0][i] for i in indices]
                    ties_counts = [ties[1][i] for i in indices]
                    if (len(ties_numbers) > 0):
                        for i in range(len(ties_numbers)):
                            #ties$count[i] * (ties$numTies[i]^3 - ties$numTies[i])
                            correction += (ties_counts[i] * (ties_numbers[i]**3 - ties_numbers[i]))
                        correction = correction / (12 * (len(tmp_values) - 1))
                weights = []
                for i in range(len(ranks[1])):
                    if (ranks[1][i] == 1):
                        # singleton, just add up previous number of values to get rank
                        weights.append(sum(ranks[1][0:i + 1]))
                    else:
                        start = sum(ranks[1][0:i])
                        ranknum = 0
                        for j in range(1, ranks[1][i] + 1, 1): # 1:ranks$count[i])
                            ranknum = ranknum + start + j
                        ranknum = ranknum / ranks[1][i]
                        weights.append(ranknum)
                ranks2 = tuple([ranks[0], ranks[1], np.array(weights)]) # ['value', 'count', 'rank']
                tmp_ranks = []
                tmp_counts = []
                for i in range(len(tmp_values)):
                    r_value = 0
                    cnt = 0
                    for j in range(len(ranks2[0])):
                        if (tmp_values[i] == ranks2[0][j]):
                            r_value = ranks2[2][j]
                            cnt = ranks[1][j]
                            break
                    tmp_ranks.append(r_value)
                    tmp_counts.append(cnt)
                tmp_data = pd.DataFrame({'doses':tmp_doses,
                                         'values':tmp_values,
                                         'ranks':tmp_ranks,
                                         'counts': tmp_counts})
                rankSums = tmp_data.groupby('doses')['ranks'].agg(['sum','count'])
                rankMeans = tmp_data.groupby('doses')['ranks'].agg(['mean','count'])
                rankSums['index'] = range(len(rankSums))
                rankSums = rankSums.reset_index().set_index('index')
                mean_line = []
                numer, denom = (0, 0)
                for i in range(len(rankSums) - 1, -1, -1):
                    numer = numer + rankSums.loc[i, 'sum']
                    denom = denom + rankSums.loc[i, 'count']
                    if (i == 0):
                        mean_temp = rankSums.loc[i, 'sum'] / rankSums.loc[i, 'count']
                    else:
                        mean_temp = numer / denom
                    mean_line.append(mean_temp)
                rankMeans['mean'] = mean_line[::-1]
                # find test statistic
                V  = (n_total * (n_total + 1)) / 12 - correction
                Ri = doses_ranks[1][-1]
                C  = doses_ranks[1][0]
                trt_means = rankMeans['mean'].tolist()
                zero_mean = trt_means.pop(0)
                if(tau >= 0):
                    dosemean = max(trt_means)
                    shrl_num = dosemean - zero_mean
                else:
                    dosemean = min(trt_means)
                    shrl_num = zero_mean - dosemean
                T = shrl_num * (V * (1/Ri + 1/C))** (-0.5)	## shirlstat in SAS code
                test_stats.append(T)
                dose_count.append(g + 1)
                test_doses.append(dose_groups[g + 1])
                test_nums.append(doses_ranks[1][g + 1])
            tests = _mult_comparison(dose_count, test_nums, test_stats, test_doses)
            results.update(tests)
        else:
            warn_message = 'Either no enough treatment groups or the control group is missing'
            results.update({'has_errors':bool(1),
                            'warnings':[warn_message]})
    return results
