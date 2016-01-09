# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from numpy import vstack, size
from time import sleep

def report_protein(protein, p, coef, i, l):
    print("Current: {} Degrad Coeff: {} Degrad: {} Coef: {} New: {} ".format(protein[i], p[0], p[0] * protein[i], coef, coef * p[l]))
    sleep(0.01)
    
# simulate the system with parameters p
def simulate(p, plt = plt, title = '', test = False):
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
        
        lacI_inh_lacI = p[9]
        if p[9] == 0:
            lacI_inh_lacI = 1 / (1 + p[4] * protein_levels[0])
        
        lacI_inh_tetR = p[9]
        if p[9] == 0:
            lacI_inh_tetR = 1 / (1 + p[5] * protein_levels[0])    
        
        tetR_inh_cI = p[8]
        if p[8] == 0:
            tetR_inh_cI = 1 / (1 + p[6] * protein_levels[1])
        
        
        # LacI
        # First term is previous protein amount
        # Second term is protein degradation depending on current amount
        # Third term is LacI inhibition on LacI default gene expression, depending on LacI amount
        protein_levels_new[0] = protein_levels[0] + p[0] * protein_levels[0] + lacI_inh_lacI * p[2];
        #TetR
        protein_levels_new[1] = protein_levels[1] + p[0] * protein_levels[1] + lacI_inh_tetR * p[2]; # calculating next protein level depending on the previous ones and parameters.
        #λcI
        protein_levels_new[2] = protein_levels[2] + p[0] * protein_levels[2] + tetR_inh_cI * p[3]; # calculating next protein level depending on the previous ones and parameters.
        #YFP
        protein_levels_new[3] = protein_levels[3] + p[0] * protein_levels[3] + (1 / (1 + p[7] * protein_levels[2])) * p[3]; # calculating next protein level depending on the previous ones and parameters.
        
        
        if test:
            protein_levels_plot = vstack((protein_levels_plot, protein_levels_new)) 
        
        protein_levels = protein_levels_new
        
        
        # if we're at a full hour, save the value
        # (that's the measured data we have and want to compare to)
        if (time == (4*60*60) or time == (6*60*60) or time == (8*60*60) or time == (10*60*60)) == 0:
            yfp_levels.append(protein_levels_new[3])
    
    if test:
        n_rows = total_time / step
        x = range(0, n_rows + 1)
        plt.plot(x, protein_levels_plot[:,0], 'r', label='LacI')
        plt.plot(x, protein_levels_plot[:,1], 'b', label='TetR')
        plt.plot(x, protein_levels_plot[:,2], 'g', label='cI')
        plt.plot(x, protein_levels_plot[:,3], 'darkgoldenrod', label='YFP')
        plt.title(title)
        plt.yscale('log')
        plt.ylabel('Protein Amount [AU]')
        plt.xlabel('time [min]')
        
        plt.legend(loc='lower left', shadow=True, fontsize='large')
        plt.show()
    return yfp_levels


# test output
# print(simulate([1,1,1]))

