# -*- coding: utf-8 -*-
import numpy as np
import itertools
import csv
import matplotlib.pyplot as plt
from simulate import simulate
from cmath import log
from operator import truediv


class ParamEvaluator(): 
    
    def __init__(self, filename):
        # string representation of the 48 types
        self.types = []
        # if the type was clonable
        self.valids = []
        # all the 6 different adjacencies are assigned a number
        self.adjdict = {'LT' : 0, 'LC' : 1, 'TC' : 2,
                        'TL' : 0, 'CL' : 1, 'CT' : 2}

        # data: 3D array with yfp expression levels for each type
        # first index: typeid
        # second index: 0, 1, 2 or 3 for wo, atc, iptg, both
        # third index: 0, 1, 2 or 3 for after 4, 6, 8 or 10 hours
        datalist = []
        
        # Number of strains in the data file (can be different than 48 for test files, expected atc file or wt etc... 
        strain_count = 0

        with open(filename, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                #print row[0]
                self.types.append(row[0])
                strain_count += 1
                rowdata = []
                valid = True
                for entry in row[1:17]:
                    #print entry
                    if entry == "":
                        valid = False
                        break
                    else:
                        rowdata.append(float(entry))
                self.valids.append(valid)
                datalist.append([rowdata[0:16:4], rowdata[1:16:4], rowdata[2:16:4], rowdata[3:16:4]])
        
        
        self.data = np.array(datalist)
        self.strain_count = strain_count
    
    @staticmethod
    def apply_ruleset(p, typeid, iptgatc):
        orient = typeid[0:3]
        arr = typeid[3:6]
        iptg = iptgatc // 2 # integer division
        atc = iptgatc % 2   # modulo

        params = np.zeros(15)
        params[10] = 1
        params[11] = 1
        
        for i in range(0, 8):
            params[i] += p[i]
        
        # if iptg is present:
        if iptg:
            params[9] += p[9]
        # if atc is present:
            if ("LT" in arr) or ("TL" in arr):
                #print(arr)
                params[10] = p[10]
        if atc:
            params[8] += p[8]
            if iptg:
                if arr[2] != 'C':
                    params[11] = p[11]

        # supercoiling penalty for genes facing each other
        if orient[0:2] == 'FR':
            params[12 + self.adjdict[arr[0:2]]] += p[12]
        if orient[1:3] == 'FR':
            params[12 + self.adjdict[arr[1:3]]] += p[12]
        
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
        
    # p: parameters that determine the strength of the effect of the rules
    # arr: order the three genes are arranged e.g. ['L', 'T', 'C']
    # ori: their orientations
    # iptg: True if IPTG is present
    # atc: True if atc is present
    # TODO: put actual rules here
    def apply_ruleset_intern(self, p, typeid, iptgatc):
        return self.apply_ruleset(p, self.types[typeid], iptgatc)


    # get_badness is the function we want to optimize
    # p is its list (array) of parameters to compute the badness for
    def get_badness(self, p, method):
        # this is the value that will accumulate the deviation from the measurements
        badness = 0.0
        # for all 48x4 possible setups
        for typeid in range(self.strain_count):
            # only simulate the cloneable arrangements
            if self.valids[typeid] == True:
                for iptgatc in range(4):
                    # get the parameters for the simulation via teh ruleset
                    params = self.apply_ruleset_intern(p, typeid, iptgatc)
                    # get the simualted yfp levels
                    yfps = np.array(simulate(params))
                    # get the actual measurements for comparison
                    measurements = self.data[typeid][iptgatc]
                    # comute the quadratic difference and add it to the badness
                    #print(yfps)
                    #print(measurements)
                    if method == 0:
                        badness += np.sum((yfps-measurements)**2)
                    elif method == 1:
                        yfps = np.maximum(yfps, np.add(np.zeros(4), 0.000001))
                        badness += np.sum(np.abs(np.log10(yfps) - np.log10(measurements)))
                    elif method == 2:
                        badness += np.sum(abs(yfps-measurements))
                        
        return badness


# func is the function we want to optimize
# init is the array of initial values of the parameters
# mins is the array of minimum values each parameter can take
# maxs is the array of maximum values each parameter can take
# method: use square differences, log ratios or absolute differences
# debug: print intermediate values
def optimize(func, init, mins, maxs, method, debug):
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
    badness = func(vals, method)
    print("%f @ temp %f: %s" % (badness, 0.21, str(vals)))
    badness_new = 0.0
    while temp > 0.0005:
        #print("Temp: %f Badness: %f" % (temp, badness_new))
        # get the new array of parameters by sampling a normal
        # distribution around the old values with standard deviation
        # proportional to temperature and range of the parameter
        vals_new = np.random.normal(vals, temp*ranges, n)
        # make sure the new values don't exceed the range
        vals_new = np.maximum(vals_new, mins)
        vals_new = np.minimum(vals_new, maxs)
        # compute the new badness
        badness_new = func(vals_new, method)
        # pick the new set of values if they are better
        if badness_new < badness:
            badness = badness_new
            vals = vals_new
            if debug:
                print("%f @ temp %f: %s" % (badness_new, temp, str(vals_new)))
        # reduce the temperature
        temp -= 0.001

    return badness, vals
    


if __name__ == "__main__":
              # init, min, max
    params = [(  0.02,     0,    0.2), # 0 Protein degradation:
              # Purcell, Oliver, and Nigel J Savery. "Temperature dependence of ssrA-tag mediated protein degradation" Jbe 6:10 (2012)
              # Halftime approx. 20% degradation in 10mins -> 2% is a good starting point 
              (   234,     0, 10000), # 1 Pλ (YFP) default expression with no interference based on max change in fluorescence readout
              # Based on fluorescence levels. Simply fluorescence levels were taken as arbitrary unit of protein amount in the cell
              (   200,     0, 10000), # 2 PLac Needs to be explored. Leakiness of the promoter is probably derived from the 1/600 repressive effect
              (   200,     0, 10000), # 3 PTet Needs to be explored. 
              (     6,     0,    40), # 4 repressive effect of LacI on LacI per 1 AU protein, will be multiplied with the protein amount
              (   5.9,     0,    40), # 5 repressive effect of LacI on TetR per 1 AU protein, will be multiplied with the protein amount
              (    50,     0,   100), # 6 repressive effect of TetR on λcI per 1 AU protein
              (   100,     0,   100), # 7 repressive effect of λcI on Pλ (YFP) to be multiplied with default expression rate
              (   0.9,   0.5,     1), # 8 inhibitory effect of aTc on tetR
              (   0.9,   0.5,     1), # 9 inhibitory effect of IPTG on pLac
              (   0.1,     0,     1), # 10 mystery inhibitory effect of IPTG and LacI-TetR neighbourship on TetR
              (     1,  0.01,     2), # 11 mystery effect when IPTG and aTc are present and C is not in the middle on TetR
              (   0.5,   0.0,     1)  # 12 mystery effect of supercoiling (two genes facing each other)
              ]
    
    #all_types = ["FFFCLT", "FFFCTL", "FFFLCT", "FFFLTC", "FFFTCL", "FFFTLC", "FRFCLT", "FRFCTL", "FRFLCT", "FRFLTC", "FRFTCL", "FRFTLC", "FFRCLT", "FFRCTL", "FFRLCT", "FFRLTC", "FFRTCL", "FFRTLC", "FRRCLT", "FRRCTL", "FRRLCT", "FRRLTC", "FRRTCL", "FRRTLC", "RRRCLT", "RRRCTL", "RRRLCT", "RRRLTC", "RRRTCL", "RRRTLC", "RRFCLT", "RRFCTL", "RRFLCT", "RRFLTC", "RRFTCL", "RRFTLC", "RFFCLT", "RFFCTL", "RFFLCT", "RFFLTC", "RFFTCL", "RFFTLC", "RFRCLT", "RFRCTL", "RFRLCT", "RFRLTC", "RFRTCL", "RFRTLC"]
    all_types = ["FFFCLT"]
    #measurement_file = 'absolute.csv';
    measurement_file = 'expected2.csv';
    #measurement_file = 'wt.csv';
    run_optization = False
    #method = 0 # quad diff
    #method = 1 # ratio-log
    method = 2 # linear diff

    if not run_optization:
        # Parameters optimized to the expected phenotype (all atc- phenotypes)
        params = [(0.083046, 0.000000, 0.200000),  # badness: 147812.035937 adjusted to: expected2.csv
                    (1233.563963, 0.000000, 10000.000000),
                    (673.542274, 0.000000, 10000.000000),
                    (1.254526, 0.000000, 10000.000000),
                    (1.245438, 0.000000, 40.000000),
                    (0.306181, 0.000000, 40.000000),
                    (51.976855, 0.000000, 100.000000),
                    (92.996234, 0.000000, 100.000000),
                    (0.946784, 0.500000, 1.000000),
                    (0.643717, 0.500000, 1.000000),
                    (0.044964, 0.000000, 1.000000),
                    (1.031398, 0.010000, 2.000000),
                    (0.966240, 0.000000, 1.000000)]
    transpose = list(zip(*params))
    init = np.array(transpose[0])
    mins = np.array(transpose[1])
    maxs = np.array(transpose[2])
    

    if run_optization:
        pe = ParamEvaluator(measurement_file)
        badness, vals = optimize(pe.get_badness, init, mins, maxs, method, debug = True)

        print("Best badness: %f" % badness)
        print("Parameters: %s" % str(vals))
        # printing optimized values in python style to include new parameters in the code 
        for i in range(len(init)):
            if i == 0:
                print("params_opt = [(%f, %f, %f),  # badness: %f adjusted to: %s" % (vals[i], mins[i], maxs[i], badness, measurement_file))
            elif (i != 0) and (i < (len(init) - 1)):
                print("\t\t\t(%f, %f, %f)," % (vals[i], mins[i], maxs[i]))
            else:
                print("\t\t\t(%f, %f, %f)]" % (vals[i], mins[i], maxs[i]))

        init = vals

    # plot graphs for the simulation
    titles = ['None', 'aTc', 'IPTG', 'Both']
    for type in all_types:
        plt.figure()
        plt.suptitle(type, fontsize=18)
        for iptgatc in range(4):
            pos = 221 + iptgatc
            plt.subplot(pos)
            applied_params = ParamEvaluator.apply_ruleset(init, type, iptgatc)
            simulate(applied_params, titles[iptgatc], test = True)
            plt.get_current_fig_manager().resize(1000, 800)
            plt.tight_layout()
        plt.show()
        #plt.savefig(type + ".png")
