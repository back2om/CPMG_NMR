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
from itertools import cycle


font = {'family' : 'normal', 'weight' : 'bold', 'size'   : 18}
matplotlib.rc('font', **font)


kex_apo = [] 
ratio_apo = []
Angle_apo = []

dR2_apo = []
with open('../../XlnB2_Slow/800MHz/Col_Res_resnum_dR2_R2o_k_dw_pa.txt') as f:
    for line in f:
        line = line.split()
        dR2_apo.append(float(line[3]))

with open('../../XlnB2_Slow/800MHz/kex_se.txt') as f:
    for line in f:
        kex_apo.append(float(line.split()[0]))

with open('../../XlnB2_Slow/800MHz/ratio_se.txt') as f:
    for line in f:
        ratio_apo.append(float(line.split()[0]))

with open('../../XlnB2_Slow/800MHz/Gradient_first_n_last_points/Angle_se.txt') as f:
    for line in f:
        Angle_apo.append(float(line))

Col_num_apo =[]
Res_name_apo =[]
R2_vals_apo = []
Std_vals_apo = []
Res_num_apo=[]

with open('../../XlnB2_Slow/800MHz/Col_Res_resnum_R2vals_std.txt') as f:
    for line in f:
         line = line.split()
         Col_num_apo.append(int(line[0]))
         Res_name_apo.append(line[1])
         Res_num_apo.append(int(line[2]))
         R2_vals_apo.append([float(i) for i in line[3:15]])
         Std_vals_apo.append([float(i) for i in line[15::]])

#########################################################################################
######## Plot 1 ##########################################################################

Indices_apo = []
Resnum_titration = [11,17,8,120,10,83, 165, 166,176, 19, 133,124, 173,48, 55, 56,171, 16, 78, 37]
Resnum_titration = [55, 56,171, 16, 78, 37]

for val in Resnum_titration:
    Indices_apo.append(Res_num_apo.index(val))



from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot(121, projection='3d')

ax1.scatter(np.array(ratio_apo), np.array(kex_apo), np.array(Angle_apo), '.', s=15,color='k')
ax1.scatter(np.array(ratio_apo)[Indices_apo], np.array(kex_apo)[Indices_apo], np.array(Angle_apo)[Indices_apo], marker='*', s=300,color='b')


myIterator = cycle([-50000, -100000, 50000, 100000])
for i in Indices_apo:
    ax1.text(ratio_apo[i] + myIterator.next(),kex_apo[i],Angle_apo[i], '%s' % (str(Res_num_apo[i])), size=12,weight='bold', color='b')
    print i, Res_num_apo[i]


plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_xlabel('$\Phi$/K')
ax1.set_ylabel('K')
ax1.set_zlabel('$\Theta$')

ax1.set_title('apo_XlnB2')
ax1.set_xlim(-10, 800000)
ax1.set_ylim(-50, 40000)
ax1.set_zlim(-2,70)
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
kex_xylo = [] 
ratio_xylo = []
Angle_xylo = []

dR2_xylo = []
with open('../../XlnB2_xylobiose_slow/800MHz/Col_Res_resnum_dR2_R2o_k_dw_pa.txt') as f:
    for line in f:
        line = line.split()
        dR2_xylo.append(float(line[3]))

with open('../../XlnB2_xylobiose_slow/800MHz/kex_se.txt') as f:
    for line in f:
        kex_xylo.append(float(line.split()[0]))

with open('../../XlnB2_xylobiose_slow/800MHz/ratio_se.txt') as f:
    for line in f:
        ratio_xylo.append(float(line.split()[0]))

with open('../../XlnB2_xylobiose_slow/800MHz/Gradient_first_n_last_points/Angle_se.txt') as f:
    for line in f:
        Angle_xylo.append(float(line))

Col_num_xylo =[]
Res_name_xylo =[]
R2_vals_xylo = []
Std_vals_xylo = []
Res_num_xylo=[]

with open('../../XlnB2_xylobiose_slow/800MHz/Col_Res_resnum_R2vals_std.txt') as f:
    for line in f:
         line = line.split()
         Col_num_xylo.append(int(line[0]))
         Res_name_xylo.append(line[1])
         Res_num_xylo.append(int(line[2]))
         R2_vals_xylo.append([float(i) for i in line[3:15]])
         Std_vals_xylo.append([float(i) for i in line[15::]])

#########################################################################################
######## Plot 1 ##########################################################################

Indices_xylo = []

for val in Resnum_titration:
    Indices_xylo.append(Res_num_xylo.index(val))

ax2 = fig.add_subplot(122, projection='3d')

ax2.scatter(np.array(ratio_xylo), np.array(kex_xylo), np.array(Angle_xylo), '.', s=15,color='k')
ax2.scatter(np.array(ratio_xylo)[Indices_xylo], np.array(kex_xylo)[Indices_xylo], np.array(Angle_xylo)[Indices_xylo], marker='*', s=300,color='b')


myIterator = cycle([-50000, -100000, 50000, 100000])
for i in Indices_xylo:
    ax2.text(ratio_xylo[i] + myIterator.next(),kex_xylo[i],Angle_xylo[i], '%s' % (str(Res_num_xylo[i])), size=12,weight='bold', color='b')
    print i, Res_num_xylo[i]


plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_xlabel('$\Phi$/K')
ax2.set_ylabel('K')
ax2.set_zlabel('$\Theta$')
ax2.set_title('XlnB2_xylobiose')
ax2.set_xlim(-10, 800000)
ax2.set_ylim(-50, 40000)
ax2.set_zlim(-2,70)

plt.show()
############################
