# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from numpy import vstack, size
from time import sleep
from math import exp

def report_protein(protein, p, coef, i, l):
    print("Current: {} Degrad Coeff: {} Degrad: {} Coef: {} New: {} ".format(protein[i], p[0], p[0] * protein[i], coef, coef * p[l]))
    sleep(0.01)
    
# simulate the system with parameters p
def simulate(p, plot = False):
    # start with the empty list
    # these will be our hourly yfp levels
    yfp_levels = []
    
    protein_levels = np.zeros(4) # Protein levels of LacI, TetR, lambdacI and YFP (respectively), which will change in each loop
    
    protein_levels_plot = np.zeros((0, 4))
    
    total_time = 10*60*60 + 1
    step = 60
    # here time increases in 60 second steps.
    # you'll probably want that to be more fine-grained
    for time in range(0, total_time, step):
        # hypothetical yfp fluorescence level
        
        #params p needs to be decided and defined (e.g. first parameter is repressive effect of LacI on TetR)
        protein_levels_new = np.zeros(4)
        
        # These terms act on the default gene expression of the corresponding genes.
        lacI_inh_lacI = p[9] + (1 - p[9]) / (1 + p[4] * protein_levels[0])
        lacI_inh_tetR = p[9] + (1 - p[9]) / (1 + p[5] * protein_levels[0])
        tetR_inh_cI   = p[8] + (1 - p[8]) / (1 + p[6] * protein_levels[1])
        cI_inh_YFP    =         1         / (1 + p[7] * protein_levels[2])
        LT_IPTG_inh_TetR = p[10] # See if this value was set. This applies to the qPCR data where TetR is drastically repressed by IPTG.
        IPTG_aTc_IFXNOR_inh_TetR = p[11] # See if this value was set. This applies to the IF and XNOR observations where IPTG can change if added with aTc and cI is not in the middle.        
        

        lacI_production = lacI_inh_lacI * p[2] * (1 + p[15])
        tetR_production = IPTG_aTc_IFXNOR_inh_TetR * LT_IPTG_inh_TetR * lacI_inh_tetR * p[2] * (1 + p[16])
        cI_production   = tetR_inh_cI * p[3] * (1 + p[17])
        YFP_production  = cI_inh_YFP * p[1]

        # supercoiling
        lacI_factor = 1 - p[12] * (1 - exp(- (lacI_production * tetR_production) / 1e6))
        tetR_factor = 1 - p[12] * (1 - exp(- (lacI_production * cI_production) / 1e6))
        cI_factor =   1 - p[12] * (1 - exp(- (tetR_production * cI_production) / 1e6))

        # First term is previous protein amount taking into account protein degradation
        # LacI # Second term is LacI inhibition on LacI default gene expression, depending on LacI amount
        protein_levels_new[0] = (1 - p[0]) * protein_levels[0] + lacI_production * lacI_factor
        #TetR # calculating next protein level depending on the previous ones and parameters.
        protein_levels_new[1] = (1 - p[0]) * protein_levels[1] + tetR_production * tetR_factor
        #λcI # calculating next protein level depending on the previous ones and parameters.
        protein_levels_new[2] = (1 - p[0]) * protein_levels[2] + cI_production * cI_factor
        #YFP # calculating next protein level depending on the previous ones and parameters.
        protein_levels_new[3] = (1 - p[0] / 3) * protein_levels[3] + YFP_production
        
        
        if plot:
            protein_levels_plot = vstack((protein_levels_plot, protein_levels_new)) 
        
        protein_levels = protein_levels_new
        
        
        # if we're at a full hour, save the value
        # (that's the measured data we have and want to compare to)
        if (time == (4*60*60) or time == (6*60*60) or time == (8*60*60) or time == (10*60*60)):
            #print(protein_levels_new[3])
            yfp_levels.append(protein_levels_new[3])
    
    if plot:
        n_rows = total_time / step
        x = range(0, int(n_rows) + 1)
        plt.plot(x, protein_levels_plot[:,0], 'r', label='LacI')
        plt.plot(x, protein_levels_plot[:,1], 'b', label='TetR')
        plt.plot(x, protein_levels_plot[:,2], 'g', label='cI')
        plt.plot(x, protein_levels_plot[:,3], 'darkgoldenrod', label='YFP', linewidth=2.0)

    #print(yfp_levels)
    return yfp_levels


# test output
# print(simulate([1,1,1]))

