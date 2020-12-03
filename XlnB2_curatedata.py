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

dw_se = []
kex_se = [] 
pa_se = []
r2o_opt_se = []
Del_R2 = []
phi_se = []
ratio_se = []
Dplus_se = []
Dmin_se = []


info_file = open('Col_Res_resnum_R2vals_std.txt', 'w')
parameter_file = open('Col_Res_resnum_dR2_R2o_k_dw_pa.txt', 'w')

workbook_xlnb2=xlrd.open_workbook("/Users/omchoudhary/Projects/CPMG_xB/CPMG_data/CPMG_XlnB2_rawdata.xlsx")
sheet_names = workbook_xlnb2.sheet_names()
workbook_xlnb2_500MHz = workbook_xlnb2.sheet_by_name(sheet_names[0])

#####################################################################
### xdata values are in 1/sec. R2 is also 1/sec. Therefore Kex obtained is 1/sec
xdatavalues = workbook_xlnb2_500MHz.col_values(0)[2:14]
xdatavalues = np.array(xdatavalues)
########################################################################################
################ Iterating across all residues
#### specify number of Cols to read only residues and not substrate
n_cols = 299
for col_num in range(1,n_cols,2):

    print col_num
    R2values=workbook_xlnb2_500MHz.col_values(col_num)[2:14]
    R2values = np.array(R2values)
    Std_values=workbook_xlnb2_500MHz.col_values(col_num+1)[2:14]

    ##### Write info file ##################################################################
    info_file.write('{:} {:}'.format(np.round(col_num), ' '))
    info_file.write('{:} {:}'.format(np.str(workbook_xlnb2_500MHz.row_values(0)[col_num])[0:3] , ' ')) 
    info_file.write('{:} {:}'.format(np.str(workbook_xlnb2_500MHz.row_values(0)[col_num])[3::] , ' '))

    for rv in R2values:
        info_file.write('{} {}'.format(rv, ' '))

    sigma = numpy.array(Std_values)
    for sv in sigma:
        info_file.write('{} {}'.format(sv, ' '))            
        
    info_file.write('\n')

    ############### Optimization ###########################################################
    # Initial guess for 3 unknowns R2o, phi, kex
    x0    = numpy.array([15, 500, 500,0.91])

    ### Objective function
    def func(x, R2o, kex, dw, pa):
        psi = kex*kex - dw*dw
        zeta = -2*dw*kex*(pa - (1-pa))
        rt = np.sqrt(psi*psi + zeta*zeta)
        Dplus = 0.5*(1 + ((psi+2*dw*dw)/rt))
        Dmin = 0.5*(-1 + ((psi+2*dw*dw)/rt))
        Nplus = (1/(1.414*x))*np.sqrt(rt + psi)
        Nmin =  (1/(1.414*x))*np.sqrt(rt - psi)
        #print (0.5/(x))
        #print np.sqrt(psi + np.sqrt(psi*psi + xi*xi))
        #print Dmin*np.cos(Nmin)
        #print Dplus*np.cosh(Nplus) - Dmin*np.cos(Nmin)
        #print Nplus
        return R2o + 0.5*(kex - (x)*(np.arccosh(Dplus*np.cosh(Nplus) - Dmin*np.cos(Nmin)))) 

    param_bounds = ([-np.inf, 0, 0.001, 0.51],[np.inf, np.inf, np.inf, 1])
    
    Opt_func = optimization.curve_fit(func, xdatavalues, R2values, x0, sigma=sigma, bounds=param_bounds, maxfev=3000)
    R2o_opt = Opt_func[0][0]
    kex_opt = Opt_func[0][1]
    dw_opt = Opt_func[0][2]
    pa_opt = Opt_func[0][3]
    phi_opt = pa_opt*(1-pa_opt)*dw_opt*dw_opt

    r2o_opt_se.append(R2o_opt)
    kex_se.append(kex_opt)
    dw_se.append(dw_opt)
    pa_se.append(pa_opt)
    phi_se.append(phi_opt)
    ratio_se.append(phi_opt/kex_opt)
    ##### Write parameter file ##################################################################


    parameter_file.write('{} {}'.format(np.round(col_num), ' '))
    parameter_file.write('{} {}'.format(np.str(workbook_xlnb2_500MHz.row_values(0)[col_num])[0:3], ' '))        
    parameter_file.write('{} {}'.format(np.str(workbook_xlnb2_500MHz.row_values(0)[col_num])[3::], ' '))
    Del_R2.append(R2values[0]-R2values[-1])
    parameter_file.write('{} {}'.format(R2values[0]-R2values[-1], ' '))
    parameter_file.write('{} {}'.format(R2o_opt, ' '))
    parameter_file.write('{} {}'.format(kex_opt, ' '))
    parameter_file.write('{} {}'.format(dw_opt, ' '))
    parameter_file.write('{} \n'.format(pa_opt))


    ########################
    psi_opt = kex_opt*kex_opt - dw_opt*dw_opt
    zeta_opt = -2*dw_opt*kex_opt*(pa_opt - (1-pa_opt))
    rt_opt = np.sqrt(psi_opt*psi_opt + zeta_opt*zeta_opt)
    Dplus_opt = 0.5*(1 + ((psi_opt+2*dw_opt*dw_opt)/rt_opt))
    Dmin_opt = 0.5*(-1 + ((psi_opt+2*dw_opt*dw_opt)/rt_opt))

    Dplus_se.append(Dplus_opt)
    Dmin_se.append(Dmin_opt)

    ######################################################################################
    ######################################################################################
    #########################################################################################
        

r2ooptfile = open('R2o_se.txt', 'w')
np.savetxt(r2ooptfile, r2o_opt_se, delimiter =" ", fmt = "%f")
r2ooptfile.close()

kexfile = open('kex_se.txt', 'w')
np.savetxt(kexfile, kex_se, delimiter =" ", fmt = "%f")
kexfile.close()

dwfile = open('dw_se.txt', 'w')
np.savetxt(dwfile, dw_se, delimiter =" ", fmt = "%f")
dwfile.close()

pafile = open('pa_se.txt', 'w')
np.savetxt(pafile, pa_se, delimiter =" ", fmt = "%f")
pafile.close()

phifile = open('phi_se.txt', 'w')
np.savetxt(phifile, phi_se, delimiter =" ", fmt = "%f")
phifile.close()

ratiofile = open('ratio_se.txt', 'w')
np.savetxt(ratiofile, ratio_se, delimiter =" ", fmt = "%f")
ratiofile.close()

Dplusfile = open('Dplus_se.txt', 'w')
np.savetxt(Dplusfile, Dplus_se, delimiter =" ", fmt = "%f")
Dplusfile.close()

Dminfile = open('Dmin_se.txt', 'w')
np.savetxt(Dminfile, Dmin_se, delimiter =" ", fmt = "%f")
Dminfile.close()


DelR2file = open('DelR2_se.txt', 'w')
np.savetxt(DelR2file, Del_R2, delimiter =" ", fmt = "%f")
DelR2file.close()

xdval = open('Xdata_Values.txt', 'w')
np.savetxt(xdval, xdatavalues, delimiter =" ", fmt = "%f")
xdval.close()


info_file.close()
parameter_file.close()
