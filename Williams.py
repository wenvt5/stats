# -*- coding: utf-8 -*-
"""
Created on Mon Aug 2 10:38:45 2018

@author: yingw
"""
import pandas as pd
import numpy as np
import os
import warnings
import williams_criticals as will



__all__ = ['Williams']


def Williams(x,y,jonck_trend):
    if (len(x) == len(y)):
        will005 = will.get_will005_csv()
        will025 = will.get_will025_csv()
        willtables = will005.merge(will025, on='dof')
        results = {'method':"Williams test",
                   'is_finished':bool(0),
                   'has_output':bool(0),
                   'has_errors':bool(0)}
        x_y_dataframe = pd.DataFrame([x,y]).transpose()
        x_y_dataframe.columns = ['x','y']
        
        ## get william-ized dose means
        ## direction of smoothing dependent on Jonckheere output
        wmeans = x_y_dataframe.groupby('x')['y'].agg(['mean','count'])
        wmeans = wmeans.reset_index(drop=False)
        wmeans = wmeans.sort_values(['x'])
        smeans = pd.Series(wmeans['mean']).copy()
        slength = pd.Series(wmeans['count']).copy()
        
        #smeans = smeans[[0,2,1,3]].reset_index(drop=True) #for test purpose
        smeans.pop(0)   ## remove control group
        slength.pop(0)  ## remove control group
        
        ## set comparison direction based on JONCK trend result
        direction = 'decreasing' if jonck_trend < 0 else 'increasing'
        
        ## need more than 1 trt grp for smoothing
        if len(smeans) > 1:
            if direction == 'decreasing':
                for i in range(1,len(smeans),1):
                    if smeans[i] < smeans[i+1]:  	## smoothing required  
                        if i==1:  ## FIRST ITEM
                            tempSmean = smeans[i]*(slength[i]/(slength[i]+slength[i+1])) + smeans[i+1]*slength[i+1]/(slength[i]+slength[i+1])
                            smeans[i] = tempSmean
                            smeans[i+1] = tempSmean
                        elif i>1: ## NOT FIRST ITEM
                            tempSmean = smeans[i]*(slength[i]/(slength[i]+slength[i+1])) + smeans[i+1]*slength[i+1]/(slength[i]+slength[i+1])
                            s_count = 0   ## initialize counter for previous consecutive smoothed means
                            for j in range(i-1, 0, -1):   	## loop back to find number of elements in smoothing
                                if tempSmean > smeans[j]:
                                    s_count = s_count + 1  ## count once for each consecutive previous "smooth"	
                                    if s_count > 0:
                                        tempSmean = 0
                                        for k in range(i-s_count, i, 1): ## recalculate tempSmean
                                            tempSmean = tempSmean + smeans[k]*slength[k]/sum(slength[i-s_count-1:i])
                                else:
                                    break
                            
                            denom = 0  	## this is denominator for weighting by sample size
                            for k in range(i+1, i-s_count-1,-1):  ## get count from previous smoothing
                                denom = denom + slength[k]
                            
                            tempSmean = 0
                            for k in range(i+1,i-s_count-1,-1): ## go back through and get smoothed mean(s)
                                tempSmean = tempSmean = smeans[k]*slength[k]/denom
                                
                            smeans[i] = tempSmean
                            smeans[i+1] = tempSmean
                            if s_count > 0:
                                for k in range(1,s_count+1,1):
                                    smeans[i-k] = tempSmean
            else:  	## DIRECTION BREAK
                for i in range(1,len(smeans),1):
                    if smeans[i] > smeans[i+1]:  ## smoothing required
                        if i==1:  	## FIRST ITEM
                            tempSmean = smeans[i]*(slength[i]/(slength[i]+slength[i+1])) + smeans[i+1]*slength[i+1]/(slength[i]+slength[i+1])
                            smeans[i] = tempSmean
                            smeans[i+1] = tempSmean
                        elif i>1: 	## NOT FIRST ITEM
                            tempSmean = smeans[i]*(slength[i]/(slength[i]+slength[i+1])) + smeans[i+1]*slength[i+1]/(slength[i]+slength[i+1])
                            s_count = 0  ## initialize counter for previous consecutive smoothed means
                            for j in range(i-1, 0, -1):  ## loop back to find number of elements in smoothing
                                if tempSmean < smeans[j]:
                                    s_count = s_count + 1  	## count once for each consecutive previous "smooth"
                                    if s_count > 0:
                                        tempSmean = 0
                                        for k in range(i-s_count, i, 1):  ## recalculate tempSmean
                                            tempSmean = tempSmean + smeans[k]*slength[k]/sum(slength[i-s_count-1:i])
                                else:
                                    break
                            
                            denom = 0  	## this is denominator for weighting by sample size
                            for k in range(i+1, i-s_count-1,-1):  ## get count from previous smoothing
                                denom = denom + slength[k]
                            
                            tempSmean = 0
                            for k in range(i+1,i-s_count-1,-1):  ## go back through and get smoothed mean(s)
                                tempSmean = tempSmean + smeans[k]*slength[k]/denom
                                
                            smeans[i] = tempSmean
                            smeans[i+1] = tempSmean
                            if s_count > 0:
                                for k in range(1,s_count+1,1):
                                    smeans[i-k] = tempSmean
                                    
        smeans = pd.concat([pd.Series(wmeans['mean'].iloc[0]),smeans]) ## pre-pend control mean to beginning of smoothed means
        wmeans['smeans'] = smeans  	## combine with means info
        
        		## get DOF for each sex/phase_type combination
			                           
		##	dof1 <- nrow(subset(will_data, sex==william$sex[w] & endpoint==william$endpoint[w] & selection==william$selection[w] & litter_name==william$litter_name[w] & phase_type==william$phase_type[w] & phase_time==william$phase_time[w] & time_in_study==william$time_in_study[w]))
		##	dof2 <- length(unique(subset(will_data, sex==william$sex[w] & endpoint==william$endpoint[w] & selection==william$selection[w] & litter_name==william$litter_name[w] & phase_type==william$phase_type[w] & phase_time==william$phase_time[w] & time_in_study==william$time_in_study[w])$dose))
		
			## simplify dof calcs ... if errors try old method above
        dof1 = sum(wmeans['count'])
        dof2 = wmeans.shape[0]
        wmeans['dof'] = dof1-dof2
        
        se = x_y_dataframe.groupby('x')['y'].agg(['std','count'])
        se['sterr'] = se['std']/(se['count']**(0.5))
        se = se.reset_index(drop=False)
        
        mse = 0
        for i in range(se.shape[0]):
            temp = (se['sterr'].iloc[i]**2) * (x_y_dataframe[x_y_dataframe['x']==se['x'].iloc[i]].shape[0]) * (x_y_dataframe[x_y_dataframe['x']==se['x'].iloc[i]].shape[0]-1)
            mse = mse + temp
            
        mse = mse/(x_y_dataframe.shape[0] - len(np.unique(x_y_dataframe['x'])))
        
        ## create WILLIAMS TEST STATISTIC
        willStat = pd.Series()
        for i in range(1,wmeans.shape[0],1):
            control_num = wmeans['count'].iloc[0]
            control_mean = wmeans['smeans'].iloc[0]
            test_num = wmeans['count'].iloc[i]
            test_mean = wmeans['smeans'].iloc[i]
            
            willstat = (control_mean - test_mean) / ((mse*((1/test_num) + (1/control_num)))**0.5)
            willStat = willStat.append(pd.Series(willstat))
        
        willStat = willStat.tolist() 
        willStat.insert(0,'control')
        wmeans['willStat'] = willStat
        
        
        	## ----------------------------------------------------------------------
		## convert williams statistic into p-value based on SAS crit levels
		## read in critical tables
		## ----------------------------------------------------------------------

		## if DOF matches with the crit tables, we can make a simple comparison, if not we extrapolate
        for k in range(wmeans.shape[0]):
            ## CONTROL GROUP
            if wmeans['x'].iloc[k] == 0:
                crit05 = 0
                crit01 = 0
                temp = pd.Series([crit05, crit01])
                
            ## NOT CONTROL AND DOF PRESENT IN TABLE
            elif wmeans['x'].iloc[k] != 0 and (wmeans['dof'].iloc[k] in willtables['dof']):
                col1 = 'w1crit' + str(k+1)
                col5 = 'w5crit' + str(k+1)
                adj1 = 'w1adj' + str(k+1)
                adj5 = 'w5adj' + str(k+1)
                
                w1crit = willtables[willtables['dof']==wmeans['dof'].iloc[k]][col1]
                w1adj = willtables[willtables['dof']==wmeans['dof'].iloc[k]][adj1]
                w5crit = willtables[willtables['dof']==wmeans['dof'].iloc[k]][col5]
                w5adj = willtables[willtables['dof']==wmeans['dof'].iloc[k]][adj5]
                
                con_num = wmeans[wmeans.x==0]['count'].iloc[0]
                trt_num = wmeans['count'].iloc[k]
                
                crit01 = w1crit - (.1 * w1adj * (1 - (trt_num / con_num)))
                crit05 = w5crit - (.1 * w5adj * (1 - (trt_num / con_num)))
                
                temp = pd.Series([crit05, crit01])
                
            	## NOT CONTROL AND DOF NOT PRESENT IN TABLE
            elif wmeans['x'].iloc[k] != 0 and not np.isin(wmeans['dof'].iloc[k],willtables['dof']).tolist():
                col1 = 'w1crit' + str(k+1)
                col5 = 'w5crit' + str(k+1)
                adj1 = 'w1adj' + str(k+1)
                adj5 = 'w5adj' + str(k+1)
                
                ## get lower bound from table
                lowdof = max(willtables[wmeans['dof'].iloc[k] > willtables['dof']]['dof'])
                
                low_w1crit = willtables[willtables['dof']==lowdof][col1].iloc[0]
                low_w1adj = willtables[willtables['dof']==lowdof][adj1].iloc[0]
                low_w5crit = willtables[willtables['dof']==lowdof][col5].iloc[0]
                low_w5adj = willtables[willtables['dof']==lowdof][adj5].iloc[0]
                
                	## get upper bound from table
                highdof = min(willtables[wmeans['dof'].iloc[k] < willtables['dof']]['dof'])
                
                high_w1crit = willtables[willtables['dof']==highdof][col1].iloc[0]
                high_w1adj = willtables[willtables['dof']==highdof][adj1].iloc[0]
                high_w5crit = willtables[willtables['dof']==highdof][col5].iloc[0]
                high_w5adj = willtables[willtables['dof']==highdof][adj5].iloc[0]
                
                dofactor = (wmeans.dof[k] - lowdof) / (highdof - lowdof)
                con_num = wmeans[wmeans.x==0]['count'].iloc[0]
                trt_num = wmeans['count'].iloc[k]
                
                crit01 = (low_w1crit - (dofactor * (low_w1crit - high_w1crit))) - (.01 * low_w1adj * (1 - (trt_num / con_num)))
                crit05 = (low_w5crit - (dofactor * (low_w5crit - high_w5crit))) - (.01 * low_w5adj * (1 - (trt_num / con_num)))
                
            will_results = wmeans.copy()
            will_results['crit01'] = crit01
            will_results['crit05'] = crit05
            
            ## remove control
            will_results = will_results[will_results.willStat != 'control']
            will_results =  will_results.reset_index()
            
            ## determine how many asterisks each row deserves
            for r in range(will_results.shape[0]):
                if will_results['willStat'].iloc[r] >= will_results['crit01'].iloc[r]:
                    will_results.loc[r,'mult_comp_signif'] = 2
                elif will_results['willStat'].iloc[r] >= will_results['crit05'].iloc[r]:
                    will_results.loc[r,'mult_comp_signif'] = 1
                else:
                    will_results.loc[r,'mult_comp_signif'] = 0
        tests = {'is_finished':bool(1),
                   'has_output':bool(1),
                   'output':{'x':will_results['x'].tolist(),
                             'means':will_results['mean'].tolist(),
                             'counts':will_results['count'].tolist(),
                             'willStats':will_results['willStat'].tolist(),
                             'mult_comp_signif':will_results['mult_comp_signif'].tolist()}}
        results.update(tests)
       
            
    else:
        warn_message = 'x and y need to be of the same length'
        results.update({'has_errors':bool(1),
                            'warnings':[warn_message]})
    return results
            
        
    
  
    
    
    
  