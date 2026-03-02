import matplotlib.pyplot as pplot
import numpy as np
import os

# USEFUL FUNCTIONS
def ReadFile(name):
    #
    return np.genfromtxt(name, delimiter=',')

# PLOT RESULTS
data = ReadFile(os.getcwd() + "/results.dat")
Figure = pplot.figure()
pplot.plot(data[:,0], data[:,1])
pplot.ylim(bottom=0)
pplot.show()