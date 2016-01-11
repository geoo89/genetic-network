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
        self.protdict = {'L' : 0, 'T' : 1, 'C' : 2}
        self.plotx = list(range(4*60, 10*60+1, 2*60))

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
    
    # p: parameters that determine the strength of the effect of the rules
    # arr: order the three genes are arranged e.g. ['L', 'T', 'C']
    # ori: their orientations
    # iptg: True if IPTG is present
    # atc: True if atc is present
    # TODO: put actual rules here
    def apply_ruleset(self, p, typeid, iptgatc):
        orient = self.types[typeid][0:3]
        arr = self.types[typeid][3:6]
        iptg = iptgatc // 2 # integer division
        atc = iptgatc % 2   # modulo

        params = np.zeros(19)
        params[10] = 1
        params[11] = 1
        
        for i in range(0, 8):
            params[i] += p[i]
        # second supercoiling param
        params[18] += p[13]
        
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
                if arr[1] != 'C':
                    params[11] = p[11]

        # supercoiling penalty for genes facing each other
        if orient[0:2] == 'FR':
            params[12 + self.adjdict[arr[0:2]]] += p[12]
        if orient[1:3] == 'FR':
            params[12 + self.adjdict[arr[1:3]]] += p[12]

        # # inverse supercoiling effect for genes facing away from each other
        # if orient[0:2] == 'RF':
        #     params[12 + self.adjdict[arr[0:2]]] += p[14]
        # if orient[1:3] == 'RF':
        #     params[12 + self.adjdict[arr[1:3]]] += p[14]

        # first gene is forward oriented -> penalty
        if orient[0] == 'F':
            params[15 + self.protdict[arr[0]]] += p[15]
        if orient[0] == 'R':
            params[15 + self.protdict[arr[0]]] += p[16]
        if orient[2] == 'F':
            params[15 + self.protdict[arr[2]]] += p[17]
        if orient[2] == 'R':
            params[15 + self.protdict[arr[2]]] += p[18]


        # at the end, return the parameters for the simulation determined by the ruleset
        return params


    # get_badness is the function we want to optimize
    # p is its list (array) of parameters to compute the badness for
    def get_badness(self, p, method, debug):
        # this is the value that will accumulate the deviation from the measurements
        badness_total = 0.0
        badnesses = []
        # for all 48x4 possible setups
        for typeid in range(self.strain_count):
            badness = 0.0
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
                    #print(yfps)
                    #print(measurements)
                    if method == 0:
                        badness += np.sum((yfps-measurements)**2)
                    elif method == 1:
                        yfps = np.maximum(yfps, np.add(np.zeros(4), 0.000001))
                        badness += np.sum(np.abs(np.log10(yfps) - np.log10(measurements)))
                    elif method == 2:
                        badness += np.sum(abs(yfps-measurements))
                    elif method == 3:
                        yfps = np.maximum(yfps, np.add(np.zeros(4), 0.000001))
                        badness += np.sum(np.exp(np.abs(np.log10(yfps) - np.log10(measurements))))
            badness_total += badness
            if debug >= 2:
                print("%s: %f" % (self.types[typeid], badness))
                if debug >= 3:
                    badnesses.append(badness)

        if debug >= 3:
            return badness_total, badnesses
        return badness_total

    def get_type(self, typeid):
        return self.types[typeid]


    def plot_measurement(self, typeid, iptgatc):
        if self.valids[typeid]:
            plt.plot(self.plotx, self.data[typeid][iptgatc], 'orange', label='YFP measured', linewidth=2.0)
    

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
    temp = 0.4
    # initial badness for the initial parameters
    badness = func(vals, method, debug)
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
        badness_new = func(vals_new, method, debug)
        # pick the new set of values if they are better
        if badness_new < badness:
            badness = badness_new
            vals = vals_new
            if debug >= 1:
                print("%f @ temp %f: %s" % (badness_new, temp, str(vals_new)))
        # reduce the temperature
        temp -= 0.002
    return badness, vals
    


if __name__ == "__main__":
              # init, min, max
    params = [(  0.02,  1e-3,   0.2), # 0 Protein degradation:
              # Purcell, Oliver, and Nigel J Savery. "Temperature dependence of ssrA-tag mediated protein degradation" Jbe 6:10 (2012)
              # Halftime approx. 20% degradation in 10mins -> 2% is a good starting point 
              (   234,    40,  3000), # 1 Pλ (YFP) default expression with no interference based on max change in fluorescence readout
              # Based on fluorescence levels. Simply fluorescence levels were taken as arbitrary unit of protein amount in the cell
              (   200,    40,  3000), # 2 PLac Needs to be explored. Leakiness of the promoter is probably derived from the 1/600 repressive effect
              (   200,    40,  3000), # 3 PTet Needs to be explored. 
              (  0.05,  1e-5,   0.1), # 4 inverse repressive effect of LacI on LacI per 1 AU protein, will be multiplied with the protein amount
              (  0.05,  1e-5,   0.1), # 5 inverse repressive effect of LacI on TetR per 1 AU protein, will be multiplied with the protein amount
              (    50,     0,   100), # 6 repressive effect of TetR on λcI per 1 AU protein
              (   100,     0,   100), # 7 repressive effect of λcI on Pλ (YFP) to be multiplied with default expression rate
              (   0.9,   0.5,     1), # 8 inhibitory effect of aTc on tetR
              (   0.9,   0.5,     1), # 9 inhibitory effect of IPTG on pLac
              (   0.1,     0,     1), # 10 mystery inhibitory effect of IPTG and LacI-TetR neighbourship on TetR
              (     1,   0.5,   1.2), # 11 mystery effect when IPTG and aTc are present and C is not in the middle on TetR
              (   0.5,   0.0,     1), # 12 mystery effect of supercoiling (two genes facing each other)
              (   1e6  , 5e5,   2e7), # 13 -> 18 mystery effect param 2 of supercoiling (two genes facing each other)
              (   0.5,  -0.2,   0.5), # 14 NOT USED mystery effect of inverse supercoiling (two genes facing away from each other)
              (   0.0,  -0.5,   0.2), # 15 mystery effect of being forward in first position (next to kanamycin resistence)
              (   0.0,  -0.5,   0.2), # 16 mystery effect of being backwd  in first position (next to kanamycin resistence)
              (   0.0,  -0.5,   0.2), # 17 mystery effect of being forward in last  position
              (   0.0,  -0.5,   0.2)] # 18 mystery effect of being backwd  in last  position
    
    #all_types = ["FFFCLT", "FFFCTL", "FFFLCT", "FFFLTC", "FFFTCL", "FFFTLC", "FRFCLT", "FRFCTL", "FRFLCT", "FRFLTC", "FRFTCL", "FRFTLC", "FFRCLT", "FFRCTL", "FFRLCT", "FFRLTC", "FFRTCL", "FFRTLC", "FRRCLT", "FRRCTL", "FRRLCT", "FRRLTC", "FRRTCL", "FRRTLC", "RRRCLT", "RRRCTL", "RRRLCT", "RRRLTC", "RRRTCL", "RRRTLC", "RRFCLT", "RRFCTL", "RRFLCT", "RRFLTC", "RRFTCL", "RRFTLC", "RFFCLT", "RFFCTL", "RFFLCT", "RFFLTC", "RFFTCL", "RFFTLC", "RFRCLT", "RFRCTL", "RFRLCT", "RFRLTC", "RFRTCL", "RFRTLC"]
    all_types = [0]
    measurement_file = 'absolute.csv';
    #measurement_file = 'expected2.csv';
    #measurement_file = 'wt.csv';
    run_optization = True
    #method = 0 # quad diff
    #method = 1 # ratio-log
    method = 2 # linear diff
    #method = 3 # exponential ratio-log

    if not run_optization:
        # Parameters optimized to the expected phenotype of FFFCTL
        params = [(0.083046, 0.000000, 0.200000),  # badness: 147812.035937 adjusted to: expected2.csv
                    (1233.563963, 0.000000, 10000.000000),
                    (673.542274, 0.000000, 10000.000000),
                    (1.254526, 0.000000, 10000.000000),
                    (6.245438, 0.000000, 40.000000),
                    (6.306181, 0.000000, 40.000000),
                    (51.976855, 0.000000, 100.000000),
                    (92.996234, 0.000000, 100.000000),
                    (0.946784, 0.500000, 1.000000),
                    (0.643717, 0.500000, 1.000000),
                    (0.044964, 0.000000, 1.000000),
                    (1.031398, 0.010000, 2.000000),
                    (0.966240, 0.000000, 1.000000)]
        

        params_opt = [(0.003831, 0.001000, 0.200000),  # badness: 6278389.546899 adjusted to: absolute.csv
            (4300.051621, 40.000000, 10000.000000),
            (6529.613055, 40.000000, 10000.000000),
            (40.000000, 40.000000, 10000.000000),
            (0.036860, 0.000010, 0.100000),
            (0.041289, 0.000010, 0.100000),
            (5.000000, 5.000000, 100.000000),
            (5.000000, 5.000000, 100.000000),
            (0.910836, 0.500000, 1.000000),
            (0.870547, 0.500000, 1.000000),
            (0.000000, 0.000000, 1.000000),
            (0.010000, 0.010000, 2.000000),
            (0.815031, 0.000000, 1.000000)]
                    
        # Parameters optimized to the expected phenotype of all types
        params = [(0.055431, 0.000000, 0.200000),  # badness: 6857086.976968 adjusted to: absolute.csv
                    (761.066718, 40.000000, 10000.000000),
                    (7739.050840, 40.000000, 10000.000000),
                    (9703.099295, 40.000000, 10000.000000),
                    (0.010490, 0.000000, 40.000000),
                    (0.010000, 0.000000, 40.000000),
                    (70.148061, 0.000000, 100.000000),
                    (22.136042, 0.000000, 100.000000),
                    (0.800482, 0.500000, 1.000000),
                    (0.760173, 0.500000, 1.000000),
                    (0.450252, 0.000000, 1.000000),
                    (0.179741, 0.010000, 2.000000),
                    (0.813581, 0.000000, 1.000000),
                    (   1e6  , 2e5,   5e7), # 13 mystery effect param 2 of supercoiling (two genes facing each other)
                    (   0.5,  -0.2,   0.5), # 14 NOT USED mystery effect of inverse supercoiling (two genes facing away from each other)
                    (   0.0,  -0.5,   0.2), # 15 mystery effect of being forward in first position (next to kanamycin resistence)
                    (   0.0,  -0.5,   0.2), # 16 mystery effect of being backwd  in first position (next to kanamycin resistence)
                    (   0.0,  -0.5,   0.2), # 17 mystery effect of being forward in last  position
                    (   0.0,  -0.5,   0.2)] # 18 mystery effect of being backwd  in last  position

        params_opt = [(0.171439, 0.001000, 0.200000),  # badness: 6775778.493931 adjusted to: absolute.csv
                    (1776.275058, 40.000000, 3000.000000), #1
                    (1609.970957, 40.000000, 3000.000000), #2
                    (40.000000, 40.000000, 3000.000000), #3
                    (0.100000, 0.000010, 0.100000), #4
                    (0.000380, 0.000010, 0.100000), #5
                    (24.456996, 0.000000, 100.000000), #6
                    (73.701940, 0.000000, 100.000000), #7
                    (0.870426, 0.500000, 1.000000), #8
                    (0.510380, 0.500000, 1.000000), #9
                    (0.033169, 0.000000, 1.000000), #10
                    (0.595020, 0.500000, 1.200000), #11
                    (0.130628, 0.000000, 1.000000), #12
                    (3982728.440083, 500000.000000, 20000000.000000), #13
                    (0.418130, -0.200000, 0.500000), #14
                    (-0.472260, -0.500000, 0.200000), #15
                    (-0.251545, -0.500000, 0.200000), #16
                    (0.025289, -0.500000, 0.200000), #17
                    (-0.225124, -0.500000, 0.200000)]


    transpose = list(zip(*params))
    init = np.array(transpose[0])
    mins = np.array(transpose[1])
    maxs = np.array(transpose[2])


    pe = ParamEvaluator(measurement_file)
    if run_optization:
        badness, vals = optimize(pe.get_badness, init, mins, maxs, method, debug = 1)

        print("Best badness: %f" % badness)
        print("Parameters: %s" % str(vals))
        # printing optimized values in python style to include new parameters in the code 
        for i in range(len(init)):
            if i == 0:
                print("params_opt = [(%f, %f, %f),  # badness: %f adjusted to: %s" % (vals[i], mins[i], maxs[i], badness, measurement_file))
            elif (i != 0) and (i < (len(init) - 1)):
                print("\t\t\t(%f, %f, %f), #%d" % (vals[i], mins[i], maxs[i], i))
            else:
                print("\t\t\t(%f, %f, %f)]" % (vals[i], mins[i], maxs[i]))

        init = vals

    # plot graphs for the simulation
    titles = ['None', 'aTc', 'IPTG', 'Both']
    bd, bdlist = pe.get_badness(init, method, debug = 3)
    for typeid in range(48):
        plt.figure(figsize=(12.0, 8.0))
        plt.suptitle("%s, badness: %.2f" % (pe.get_type(typeid), bdlist[typeid]), fontsize=18)
        for iptgatc in range(4):
            pos = 221 + iptgatc
            plt.subplot(pos)
            applied_params = pe.apply_ruleset(init, typeid, iptgatc)
            simulate(applied_params, plot = True)
            pe.plot_measurement(typeid, iptgatc)
            plt.title(titles[iptgatc])
            plt.yscale('log')
            plt.ylim([1, 100000])
            plt.ylabel('Protein Amount [AU]')
            plt.xlabel('time [min]')
            plt.legend(loc='lower right', shadow=True, fontsize='large')
        plt.get_current_fig_manager().resize(1000, 800)
        plt.tight_layout()
        #plt.show()
        plt.savefig("figures/%s.png" % pe.get_type(typeid))
        plt.close()
