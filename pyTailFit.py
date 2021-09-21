# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 10:59:33 2021

@author: Danny Petschke
@email:  danny.petschke@uni-wuerzburg.de

"""

#*************************************************************************************************
#**")
#** Copyright (c) 2021 Danny Petschke. All rights reserved.
#**")
#** This program is free software: you can redistribute it and/or modify
#** it under the terms of the GNU General Public License as published by
#** the Free Software Foundation, either version 3 of the License, or
#** (at your option) any later version.
#**
#** This program is distributed in the hope that it will be useful,
#** but WITHOUT ANY WARRANTY; without even the implied warranty of
#** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#** GNU General Public License for more details.
#**")
#** You should have received a copy of the GNU General Public License
#** along with this program. If not, see http://www.gnu.org/licenses/.
#**
#** Contact: danny.petschke@uni-wuerzburg.de
#**
#*************************************************************************************************

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func_1(x,a,tau,b): # template function for mono-exp-decay model
     return a*np.exp(-x/tau)+b
 
def func_2(x,a1,tau1,a2,tau2,b): # template function for bi-exp-decay model
     return a1*np.exp(-x/tau1)+a2*np.exp(-x/tau2)+b
 
# add additional components/models like ... 
     
# def func_3(x,a1,tau1,a2,tau2,a3,tau3,..,b):
#     return a1*np.exp(-x/tau1)+a2*np.exp(-x/tau2)+a3*np.exp(-x/tau3)+ ... +b
 
def tail_fit(spectrum=[],
             no_of_expected_decays=1,
             no_chn_right_of_peak=200,
             bin_width_in_ps=5.,
             debug=True):
    fit_start_index = int(np.argmax(spectrum)+no_chn_right_of_peak)
    time = np.linspace(0.,len(spectrum),len(spectrum))*bin_width_in_ps
    
    # apply error weightening ...
    
    fit_weighting = np.ones(len(spectrum[fit_start_index:]))
    
    for i in range(len(spectrum[fit_start_index:])):
        val = 1./np.sqrt(spectrum[fit_start_index+i]) 
        
        if np.isfinite(val): 
            fit_weighting[i] += val
        
    # fit data ...
    
    if no_of_expected_decays == 1: # mono-exp-decay ...
        p0=[0.5*np.amax(spectrum),1000.,1.] # estimate start-values
        
        popt,pcov = curve_fit(func_1,time[fit_start_index:],spectrum[fit_start_index:],maxfev=1000000,p0=p0,sigma=fit_weighting,absolute_sigma=True,bounds=(0,[1e10,142000.,1e10]))
        fit_result = func_1(time[fit_start_index:],popt[0],popt[1],popt[2]) # resulting fit curve
        
        puncertainties = np.sqrt(np.diag(pcov))
        
        if debug:
            plt.semilogy(time[fit_start_index:],spectrum[fit_start_index:],'bo',time[fit_start_index:],fit_result,'r-')
            plt.show()
        
            print('tau-1 = ({} +/- {}) ps'.format(popt[1],puncertainties[1]))
            # print('amp-1 = ({} +/- {}) #' .format(popt[0],puncertainties[0]))
            print('bkrd  = ({} +/- {}) #' .format(popt[2],puncertainties[2]))
        
        return time[fit_start_index:],spectrum[fit_start_index:],fit_result # return fit results ...
    
    elif no_of_expected_decays == 2: # bi-exp-decay ...
        p0=[0.5*np.amax(spectrum),1000.,0.5*np.amax(spectrum),1000.,1.] # estimate start-values
        
        popt,pcov = curve_fit(func_2,time[fit_start_index:],spectrum[fit_start_index:],maxfev=1000000,p0=p0,sigma=fit_weighting,absolute_sigma=True,bounds=(0,[1e10,142000.,1e10,142000.,1e10]))
        fit_result = func_2(time[fit_start_index:],popt[0],popt[1],popt[2],popt[3],popt[4]) # resulting fit curve ... 
        
        puncertainties = np.sqrt(np.diag(pcov))
        
        if debug:
            plt.semilogy(time[fit_start_index:],spectrum[fit_start_index:],'bo',time[fit_start_index:],fit_result,'r-')
            plt.show()
        
        exp_1 = func_1(time[fit_start_index:],popt[0],popt[1],0.)
        exp_2 = func_1(time[fit_start_index:],popt[2],popt[3],0.)
        
        I_1 = sum(exp_1)/(sum(exp_1)+sum(exp_2))
        I_2 = sum(exp_2)/(sum(exp_1)+sum(exp_2))
        
        if debug:
            print('tau-1 = ({} +/- {}) ps'.format(popt[1],puncertainties[1]))
            print('amp-1 = ({} +/- {}) # [>> {} %]' .format(popt[0],puncertainties[0],100.*I_1))
            print('tau-2 = ({} +/- {}) ps'.format(popt[3],puncertainties[3]))
            print('amp-2 = ({} +/- {}) # [>> {} %]' .format(popt[2],puncertainties[2],100.*I_2))
            print('bkrd  = ({} +/- {}) #' .format(popt[4],puncertainties[4]))
        
        return time[fit_start_index:],spectrum[fit_start_index:],fit_result # return fit results ...