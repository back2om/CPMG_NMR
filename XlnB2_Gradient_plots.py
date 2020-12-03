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

######### READ Data Files #########################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('residues', nargs='+', type = int,help='list of residues')
args = parser.parse_args()


xdata =[]
with open('../Xdata_Values.txt') as f:
    for line in f:
        xdata.append(float(line.strip()))
        

Del_R2 =[]
phi_se =[]
R2o_se =[]
phi_se = []
kex_se = [] 
dw_se = []
pa_se = []
ratio_se = []


with open('../kex_se.txt') as f:
    for line in f:
        kex_se.append(float(line.split()[0]))

with open('../phi_se.txt') as f:
    for line in f:
        phi_se.append(float(line.split()[0]))

with open('../R2o_se.txt') as f:
    for line in f:
        R2o_se.append(float(line.split()[0]))

with open('../ratio_se.txt') as f:
    for line in f:
        ratio_se.append(float(line.split()[0]))

with open('../Col_Res_resnum_dR2_R2o_k_dw_pa.txt') as f:
    for line in f:
        line = line.split()
        Del_R2.append(float(line[3]))
        dw_se.append(float(line[6]))
        pa_se.append(float(line[7]))


Col_num =[]
Res_name =[]
Res_num=[]
R2_vals = []
Std_vals = []

with open('../Col_Res_resnum_R2vals_std.txt') as f:
    for line in f:
        line = line.split()
        Col_num.append(int(line[0]))
        Res_name.append(line[1])
        Res_num.append(int(line[2]))
        R2_vals.append([float(i) for i in line[3:15]])
        Std_vals.append([float(i) for i in line[15::]])
 


G_lines = [line.split()[2:12] for line in open('Gradients_se.txt')]
A_lines = [line.split()[2:12] for line in open('Angles_se.txt')]

#########################################################################################
fig, axs = plt.subplots(2,4,sharex=True)
fig.subplots_adjust(bottom=0.2)
axs = axs.ravel()

#### Plot resdidues with ratio > 100000 in ARM 1 ######
plot_residues = args.residues
Indices =  [Res_num.index(i) for i in plot_residues]

print plot_residues 
print Indices
x_grid = list(arange(100, 1700,200)) 

def func(x, R2o, kex, dw, pa):
        psi = kex*kex - dw*dw
        zeta = -2*dw*kex*(pa - (1-pa))
        rt = np.sqrt(psi*psi + zeta*zeta)
        Dplus = 0.5*(1 + ((psi+2*dw*dw)/rt))
        Dmin = 0.5*(-1 + ((psi+2*dw*dw)/rt))
        Nplus = (1/(1.414*x))*np.sqrt(rt + psi)
        Nmin =  (1/(1.414*x))*np.sqrt(rt - psi)
        return R2o + 0.5*(kex - (x)*(np.arccosh(Dplus*np.cosh(Nplus) - Dmin*np.cos(Nmin)))) 


for count,k in enumerate(Indices):
    Residue_num = Res_num[k]
    Residue_name = Res_name[k]
    
    R2values_optimized = []
    for xv in x_grid:
        R2values_optimized.append(func(xv, R2o_se[k], kex_se[k], dw_se[k], pa_se[k]))

    ######## Plot 1 ##########################################################################
    ax=axs[count]
    print count
    ax.plot(x_grid, R2values_optimized, 'k', label='$pa$ {:.2f} \n$K$ {:.2f} \n$dw$ {:.2f} \n$\Delta r2$ {:.2f} \n ratio {:.2f}'.format(pa_se[k], kex_se[k], dw_se[k], Del_R2[k],ratio_se[k]))
    ax.errorbar(xdata, R2_vals[k], Std_vals[k],fmt='o',capsize=2, elinewidth=2, markeredgewidth=2, color='g')
    ax.text(.1,.9,'{}{}'.format(Residue_name, Residue_num),horizontalalignment='left', transform=ax.transAxes)
    ax.legend( prop={'size': 12},frameon=False)
    ax.grid()
    ax.set_ylim(6, 26)

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

bar_xcoor = [100, 200, 300, 400, 500, 600, 800, 1000, 1400, 1600]
for count,k in enumerate(Indices):
    ######## Plot 2 ##########################################################################
    ax=axs[count+4]
    print count
    ax.bar(bar_xcoor, [float(i) for i in G_lines[k]], align = 'center', width=100)
    ax.grid()
    ax.set_ylim(-20, 10)
    ax.set_ylabel('Gradient')
    ax.set_xticklabels([])


plt.show()
############################













