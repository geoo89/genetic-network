import numpy as np

# simulate the system with parameters p
def simulate(p):
    # start with the empty list
    # these will be our hourly yfp levels
    yfp_levels = []

    # here time increases in 60 second steps.
    # you'll probably want that to be more fine-grained
    for time in range(0, 24*60*60, 60):
        # hypthetical yfp fluorescence level
        # here you'd have to apply a simulation step
        yfp = np.sqrt(float(time)/60)*p[0] - p[1] + p[2]*(p[2]+1)
        # if we're at a full hour, save the value
        # (that's the measured data we have and want to compare to)
        if (time % (60*60)) == 0:
            yfp_levels.append(yfp)

    return yfp_levels


# test output
# print(simulate([1,1,1]))

