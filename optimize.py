import numpy as np
import itertools
from simulate import simulate


# p: parameters that determine the strength of the effect of the rules
# arr: order the three genes are arranged
# ori: their orientations
# iptg: True if IPTG is present
# atc: True if atc is present
# TODO: put actual rules here
def apply_ruleset(p, arr, ori, iptg, atc):
    params = np.zeros(3)
    # if lacI is first:
    if arr[0] == 'L':
        params[0] += p[0]
    # if lacI is right before tetR:
    if (arr[0] == 'L' and arr[1] == 'T') or (arr[1] == 'L' and arr[2] == 'T'):
        params[1] += p[1]

    # at the end, return the parameters for the simulation determined by the rulset
    return params


# f is the function we want to optimize
# p is its list (array) of parameters
def f(p):
    # this is the value that will accumulate the deviation from the measurements
    badness = 0.0
    # for all 48x4 possible setups
    # TODO: optimize this into a single loop with a lookup
    for arrangement in itertools.permutations(['L', 'T', 'C']):
        for first in ['F', 'R']:
            for second in ['F', 'R']:
                for third in ['F', 'R']:
                    for iptg in [True, False]:
                        for atc in [True, False]:
                            orientations = [first, second, third]
                            params = apply_ruleset(p, arrangement, orientations, iptg, atc)
                            if params != None:
                                yfps = simulate(params)
                                # TODO: compute the difference between measurement and simulation here
                                diff = sum(yfps)
                                # add it to the badness
                                badness += diff
                            else:
                                # we've figured that some are non-cloneable
                                pass

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
    params = [(  0,  -1,    1),
              (  0,  -1,    1),
              (  0,  -1,    1)]

    transpose = list(zip(*params))
    init = np.array(transpose[0])
    mins = np.array(transpose[1])
    maxs = np.array(transpose[2])

    optimize(f, init, mins, maxs)


    