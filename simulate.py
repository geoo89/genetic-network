import numpy as np

# simulate the system with parameters p
def simulate(p):
    # start with the empty list
    # these will be our hourly yfp levels
    yfp_levels = []
    
    protein_levels = np.zeros(3) # Protein levels of LacI, TetR and λcI (respectively), which will change in each loop
    
    # here time increases in 60 second steps.
    # you'll probably want that to be more fine-grained
    for time in range(0, 24*60*60, 60):
        # hypothetical yfp fluorescence level
        
        #params p needs to be decided and defined (e.g. first parameter is repressive effect of LacI on TetR)
        protein_levels_new = np.zeros(3)
        #LacI
        protein_levels_new[0] = protein_levels[0] + p[0] * protein_levels[0]; # calculating next protein level depending on the previous ones and parameters.
        #TetR
        protein_levels_new[0] = protein_levels[0] + p[0] * protein_levels[0]; # calculating next protein level depending on the previous ones and parameters.
        #λcI
        protein_levels_new[0] = protein_levels[0] + p[0] * protein_levels[0]; # calculating next protein level depending on the previous ones and parameters.
        #YFP
        
        #TODO Edit previous line and implement for 1 (TetR) and 2 (λcI)
        
        
        
        yfp = np.sqrt(float(time)/60)*p[0] - p[1] + p[2]*(p[2]+1)
        # if we're at a full hour, save the value
        # (that's the measured data we have and want to compare to)
        if (time % (60*60)) == 0:
            yfp_levels.append(yfp)

    return yfp_levels


# test output
# print(simulate([1,1,1]))

