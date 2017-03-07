import numpy as np
import pylab


data = np.loadtxt("data.out")
spline = np.loadtxt("test.out")

pylab.plot(data[:,0], data[:,1], 'r+')
pylab.plot(spline[:,0], spline[:,1])

pylab.xlim([0,10])
pylab.ylim([-5,60])

pylab.show()