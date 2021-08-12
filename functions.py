#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in Jul 2021

@author: Naomi Rankin

"""
#==============================================================================
# libraries
#==============================================================================
import os
import copy
import numpy as np
import pandas as pd
import datetime
from pathlib import Path
import shutil
import scipy.stats as sst
import scipy.special as ssp
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.cm as cm
import matplotlib.colors as mco
import matplotlib.patches as mpatches
import matplotlib.colors as mco
import matplotlib.ticker as ticker

#==============================================================================
# mutation https://www.nature.com/articles/s41467-020-19818-2 - clock rate
#==============================================================================

def agent_mut_num(rate=0.419):
#adds mutationfor every 3 transmissions
#returns: number of mutations
    mymut = np.random.poisson(rate)
    return mymut


def spot_mutation(rate, spot_type):
#one mutation for every 2-3 transmissions
#returns: new nucleotide for each spot
    conversion = ['A', 'C', 'T', 'G']
    conversion = conversion[conversion != spot_type]
    new = random.choice(conversion)
    return new

#==============================================================================
# infection 
#==============================================================================

def num_infected(inner, outer, Sus_I, Sus_J):
#inner beta, intra beta, arrays of indices for inside and outside group
#returns: arrays of index numbers for those selected
    I_inf = []
    J_inf = []
    num_in = inner*len(Sus_I)
    num_out = inner*len(Sus_J)
    for i in range(num_in):
        I_inf.append(random.choice(Sus_I))
    for o in range(num_out):
        I_inf.append(random.choice(Sus_J)
    return I_inf, J_inf
                     
#==============================================================================
# recovery 
#==============================================================================

def recovery_odd(k=1,theta=10):
#gamma functoin - time until next event
                     #k is number of events (one)
                     #theta is avg num days to infection (10)
#returns: time in days until recovery. save and then check if infection day is higher than that
    day = round(np.random.gamma(k, theta, 1)[0])
    return day