# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 22:24:45 2021

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

import matplotlib.pyplot as plt
import numpy as np

import pyTailFit as ptf

# example code ...

if __name__ == '__main__':
    __,spectrum  = np.loadtxt('.../test-data.dat', delimiter='\t', skiprows=5, unpack=True, dtype='float')
        
    plt.semilogy(spectrum[:],'bo')
    plt.show()
    
    time,data,fit=ptf.tail_fit(spectrum=spectrum[:],no_of_expected_decays=1,no_chn_right_of_peak=400,bin_width_in_ps=8.)