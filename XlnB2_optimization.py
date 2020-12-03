import matplotlib.pyplot as plt
import matplotlib
import xlrd
from xlrd.sheet import ctype_text   
import numpy as np
import mixture, math
import matplotlib.mlab as mlab
from pomegranate import *
import scipy.optimize
import scipy.stats
from numpy import *
import scipy.optimize as optimization
import matplotlib.animation as animation

font = {'family' : 'normal', 'weight' : 'bold', 'size'   : 18}
matplotlib.rc('font', **font)


R2_Trp20  = np.array([22.34844, 21.78013, 21.02926, 18.3895, 18.94025, 18.68749, 17.82647, 17.20728, 16.644, 16.72424, 17.44654, 17.44564])
xdatavalues = np.array([ 100.0,  200.0,  300.030003,  399.99984, 400.0,  500.0, 599.880024,  800.0, 1000.0, 1400.55826252, 1400.56022409, 1600.0])
#xdatavalues = xdatavalues/1000.0
sigma = np.array([0.91855, 0.7573, 0.97549, 0.94508, 0.81928, 1.01058, 1.01058, 0.86541, 0.76841, 1.15563, 0.84827, 0.89529])

kex_range = np.arange(0.01,3000,10)
phi_range = np.arange(0.01,30000,100)
R2o_range = np.arange(-2.0, 25.0, 1)

parameter_file = open('Trp20_R2o_phi_k_error.txt', 'w')

#####################################################################
############### Optimization ###########################################################
### Objective function
def func(x, R2o, phiex, kex):
    return R2o + (phiex/kex)*(1-(2*np.tanh((1.0/x)*kex*0.5)/((1.0/x)*kex)))

param_bounds = ([-np.inf, 0, 0],[np.inf, np.inf, np.inf])


for l in R2o_range:
    print l
    for m in phi_range:
        for n in kex_range:
            x0    = numpy.array([l,m,n])
            Opt_func = optimization.curve_fit(func, xdatavalues, R2_Trp20, x0, sigma=sigma, bounds=param_bounds, maxfev=3000)
            R2o_opt = Opt_func[0][0]
            phiex_opt = Opt_func[0][1]
            kex_opt = Opt_func[0][2]

            #### Compute Error    
            def func_optimized(x, R2o_opt, phiex_opt, kex_opt):
                return R2o_opt + (phiex_opt/kex_opt)*(1-(2*np.tanh((1.0/x)*kex_opt*0.5)/((1.0/x)*kex_opt)))

            R2_optimized = []
            for k in xdatavalues:
                R2_optimized.append(func_optimized(k, R2o_opt, phiex_opt, kex_opt))

            R2_optimized=np.array(R2_optimized)
            #print R2_optimized
            #### Get Goodness of fit
            Gof = np.sqrt(((R2_Trp20-R2_optimized)**2).sum())
    
            ##### Write parameter file ##################################################################
            parameter_file.write('{} {}'.format(l, ' '))
            parameter_file.write('{} {}'.format(m, ' '))
            parameter_file.write('{} {}'.format(n, ' '))
            parameter_file.write('{} \n'.format(Gof))

######################################################################################
######################################################################################
        

parameter_file.close()
