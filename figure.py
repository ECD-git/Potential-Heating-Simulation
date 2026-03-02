import matplotlib.pyplot as pplot
import numpy as np
import os

def ReadFile(name):
    return np.genfromtxt(name, delimiter=',')


data = ReadFile(os.getcwd() + "/results.dat")
Figure = pplot.figure()
pplot.plot(data[:,0], data[:,1])
pplot.show()