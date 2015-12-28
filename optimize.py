import numpy as np


# f is the function we want to optimize
# p is its list (array) of parameters
def f(p):
    return p[0] - p[1] + p[2]*(p[2]+1)


# we want to optimize the function func
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
    while temp > 0.0005:
        # get the new array of parameters by sampling a normal
        # distribution around the old values with standard deviation
        # proportional to temperature and range of the parameter
        vals_new = np.random.normal(vals, temp*ranges, n)
        # make sure the new values don't exceed the range
        vals_new = np.maximum(vals_new, mins)
        vals_new = np.minimum(vals_new, maxs)
        # pick the new set of values if they are better
        if func(vals_new) < func(vals):
            print("%f @ temp %f: %s" % (f(vals_new), temp, str(vals_new)))
            vals = vals_new
        # reduce the temperature
        temp -= 0.0001

    print("Best function value: %f" % f(vals))
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


    