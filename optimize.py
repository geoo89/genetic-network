# -*- coding: utf-8 -*-
import numpy as np
import itertools
import csv
import matplotlib.pyplot as plt
from simulate import simulate


class ParamEvaluator(): 
    
    def __init__(self, filename):
        # string representation of the 48 types
        self.types = []
        # if the type was clonable
        self.valids = []
        # data: 3D array with yfp expression levels for each type
        # first index: typeid
        # second index: 0, 1, 2 or 3 for wo, atc, iptg, both
        # third index: 0, 1, 2 or 3 for after 4, 6, 8 or 10 hours
        datalist = []

        with open(filename, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                self.types.append(row[0])
                rowdata = []
                valid = True
                for entry in row[1:17]:
                    if entry == "":
                        valid = False
                        break
                    else:
                        rowdata.append(float(entry))
                self.valids.append(valid)
                datalist.append([rowdata[0:16:4], rowdata[1:16:4], rowdata[2:16:4], rowdata[3:16:4]])

        self.data = np.array(datalist)


    # p: parameters that determine the strength of the effect of the rules
    # arr: order the three genes are arranged e.g. ['L', 'T', 'C']
    # ori: their orientations
    # iptg: True if IPTG is present
    # atc: True if atc is present
    # TODO: put actual rules here
    def apply_ruleset(self, p, typeid, iptgatc):
        arr = self.types[typeid]
        iptg = iptgatc // 2 # integer division
        atc = iptgatc % 2   # modulo

        params = np.zeros(10)
        
        for i in range(0, 7):
            params[i] += p[i]
        
        # if iptg is present:
        if iptg:
            params[9] += p[9]
        # if atc is present:
        if atc:
            params[8] += p[8]


        # mystery stuff
        # if lacI is first:
        #if arr[3] == 'L':
        #    params[2] += p[2]
        # if lacI is right before tetR:
        #if (arr[3:4] == 'LT') or (arr[4:6] == 'LT'):
        #    params[3] += p[3]
        #    params[2] += p[4]

        # at the end, return the parameters for the simulation determined by the ruleset
        return params


    # get_badness is the function we want to optimize
    # p is its list (array) of parameters to compute the badness for
    def get_badness(self, p):
        # this is the value that will accumulate the deviation from the measurements
        badness = 0.0
        # for all 48x4 possible setups
        for typeid in range(48):
            # only simulate the cloneable arrangements
            if self.valids[typeid] == True:
                for iptgatc in range(4):
                    # get the parameters for the simulation via teh ruleset
                    params = self.apply_ruleset(p, typeid, iptgatc)
                    # get the simualted yfp levels
                    yfps = np.array(simulate(params))
                    # get the actual measurements for comparison
                    measurements = self.data[typeid][iptgatc]
                    # comute the quadratic difference and add it to the badness
                    print(yfps)
                    print(measurements)
                    badness += np.sum((yfps-measurements)**2)
        return badness


# func is the function we want to optimize
# init is the array of initial values of the parameters
# mins is the array of minimum values each parameter can take
# maxs is the array of maximum values each parameter can take
def optimize(func, init, mins, maxs):
    # number of parameters
    n = len(init)
    # the range of each of the 
    ranges = maxs - mins #np.array([maxs[i] - mins[i] for i in range(n)])
    # vals will be the current set of parameter values
    vals = init
    # temp is the heat/temperature, determining how far we may vary
    # when picking the next set of parameters to try
    temp = 0.2
    # initial badness for the initial parameters
    badness = func(vals)
    print("%f @ temp %f: %s" % (badness, 0.21, str(vals)))

    while temp > 0.0005:
        print(temp)
        # get the new array of parameters by sampling a normal
        # distribution around the old values with standard deviation
        # proportional to temperature and range of the parameter
        vals_new = np.random.normal(vals, temp*ranges, n)
        # make sure the new values don't exceed the range
        vals_new = np.maximum(vals_new, mins)
        vals_new = np.minimum(vals_new, maxs)
        # compute the new badness
        badness_new = func(vals_new)
        # pick the new set of values if they are better
        if badness_new < badness:
            print("%f @ temp %f: %s" % (badness_new, temp, str(vals_new)))
            badness = badness_new
            vals = vals_new
        # reduce the temperature
        temp -= 0.01

    print("Best function value: %f" % badness)
    print("Parameters: %s" % str(vals))


if __name__ == "__main__":
              # init, min, max
    params = [( -0.02,  -0.2,     0), # 0 Protein degradation:
              # Purcell, Oliver, and Nigel J Savery. "Temperature dependence of ssrA-tag mediated protein degradation" Jbe 6:10 (2012)
              # Halftime approx. 20% degradation in 10mins -> 2% is a good starting point 
              (   234,     0, 10000), # 1 Pλ (YFP) default expression with no interference based on max change in fluorescence readout
              # Based on fluorescence levels. Simply fluorescence levels were taken as arbitrary unit of protein amount in the cell
              (   200,     0, 10000), # 2 PLac Needs to be explored. Leakiness of the promoter is probably derived from the 1/600 repressive effect
              (   200,     0, 10000), # 3 PTet Needs to be explored. 
              (     6,     0,   100), # 4 repressive effect of LacI on LacI per 1 AU protein, will be multiplied with the protein amount
              (   5.9,     0,   100), # 5 repressive effect of LacI on TetR per 1 AU protein, will be multiplied with the protein amount
              (    50,     0,   100), # 6 repressive effect of TetR on λcI per 1 AU protein
              (    40,     0,   100), # 7 repressive effect of λcI on Pλ (YFP) to be multiplied with default expression rate
              (   0.6,     0,     1), # 8 inhibitory effect of aTc on tetR
              (   0.5,     0,     1)  # 9 inhibitory effect of IPTG on LacI
              
              ]
    
    test = False
    
    transpose = list(zip(*params))
    init = np.array(transpose[0])
    mins = np.array(transpose[1])
    maxs = np.array(transpose[2])
    
    # Testing
    if test:
        plt.ion()
        plt.subplot(221)    
        init[8] = 0
        init[9] = 0
        simulate(init, plt, "None", test = test)
        plt.subplot(222)
        init[8] = 0.2
        init[9] = 0
        simulate(init, plt, "aTc", test = test)
        plt.subplot(223)
        init[8] = 0
        init[9] = 0.1
        simulate(init, plt, "IPTG", test = test)
        plt.subplot(224)
        init[8] = 0.2
        init[9] = 0.1
        simulate(init, plt, "Both", test = test)
        plt.get_current_fig_manager().resize(1000, 800)
        plt.tight_layout()
        sleep(10)
    else:
        pe = ParamEvaluator('absolute.csv')
        optimize(pe.get_badness, init, mins, maxs)
    

    