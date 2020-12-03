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
import matplotlib.pyplot as plt
fig = plt.figure()

######### READ Data Files #########################
xdata_PpR3 =[]
with open('../PpR3/Fast_Exchange/Xdata_Values.txt') as f:
    for line in f:
        xdata_PpR3.append(float(line.strip()))
        
R2o_opt_PpR3 =[]
with open('../PpR3/Fast_Exchange/R2o_fe.txt') as f:
    for line in f:
        R2o_opt_PpR3.append(float(line.strip()))

Phiex_opt_PpR3 =[]
with open('../PpR3/Fast_Exchange/phiex_fe.txt') as f:
    for line in f:
        Phiex_opt_PpR3.append(float(line.strip()))

kex_opt_PpR3 =[]
with open('../PpR3/Fast_Exchange/Kex_fe.txt') as f:
    for line in f:
        kex_opt_PpR3.append(float(line.strip()))

Del_R2_PpR3 =[]
with open('../PpR3/Fast_Exchange/DelR2_fe.txt') as f:
    for line in f:
        Del_R2_PpR3.append(float(line.strip()))

ratio_phi_kex_PpR3 =[]
with open('../PpR3/Fast_Exchange/ratio_phiex_kex_fe.txt') as f:
    for line in f:
        ratio_phi_kex_PpR3.append(float(line.strip()))

#############################

xdata_PaR3 =[]
with open('../PaR3/Fast_Exchange/Xdata_Values.txt') as f:
    for line in f:
        xdata_PaR3.append(float(line.strip()))
        
R2o_opt_PaR3 =[]
with open('../PaR3/Fast_Exchange/R2o_fe.txt') as f:
    for line in f:
        R2o_opt_PaR3.append(float(line.strip()))

Phiex_opt_PaR3 =[]
with open('../PaR3/Fast_Exchange/phiex_fe.txt') as f:
    for line in f:
        Phiex_opt_PaR3.append(float(line.strip()))

kex_opt_PaR3 =[]
with open('../PaR3/Fast_Exchange/Kex_fe.txt') as f:
    for line in f:
        kex_opt_PaR3.append(float(line.strip()))

Del_R2_PaR3 =[]
with open('../PaR3/Fast_Exchange/DelR2_fe.txt') as f:
    for line in f:
        Del_R2_PaR3.append(float(line.strip()))

ratio_phi_kex_PaR3 =[]
with open('../PaR3/Fast_Exchange/ratio_phiex_kex_fe.txt') as f:
    for line in f:
        ratio_phi_kex_PaR3.append(float(line.strip()))

#################################

xdata_MfR3 =[]
with open('../MfR3/Fast_Exchange/Xdata_Values.txt') as f:
    for line in f:
        xdata_MfR3.append(float(line.strip()))
        
R2o_opt_MfR3 =[]
with open('../MfR3/Fast_Exchange/R2o_fe.txt') as f:
    for line in f:
        R2o_opt_MfR3.append(float(line.strip()))

Phiex_opt_MfR3 =[]
with open('../MfR3/Fast_Exchange/phiex_fe.txt') as f:
    for line in f:
        Phiex_opt_MfR3.append(float(line.strip()))

kex_opt_MfR3 =[]
with open('../MfR3/Fast_Exchange/Kex_fe.txt') as f:
    for line in f:
        kex_opt_MfR3.append(float(line.strip()))

Del_R2_MfR3 =[]
with open('../MfR3/Fast_Exchange/DelR2_fe.txt') as f:
    for line in f:
        Del_R2_MfR3.append(float(line.strip()))

ratio_phi_kex_MfR3 =[]
with open('../MfR3/Fast_Exchange/ratio_phiex_kex_fe.txt') as f:
    for line in f:
        ratio_phi_kex_MfR3.append(float(line.strip()))

##########################################################################################
######## Plot 1 ##########################################################################
ax1=fig.add_subplot(341)
ax1.scatter(kex_opt_PpR3, ratio_phi_kex_PpR3, marker='o', s=40, c='b', label = 'PpR3')
ax1.set_ylabel('$\Phi/K_{ex}$')
ax1.set_ylim([0,40])
ax1.set_xlim([-2,30])
ax1.grid()
plt.legend(prop={'size': 14})


ax2=fig.add_subplot(342)
ax2.scatter(Phiex_opt_PpR3, ratio_phi_kex_PpR3, marker='o', s=40, c='b', label = 'PpR3')
ax2.set_ylim([0,40])
ax2.set_xlim([-10,300])
ax2.grid()
plt.legend(prop={'size': 14})

ax3=fig.add_subplot(343)
ax3.scatter(R2o_opt_PpR3, ratio_phi_kex_PpR3, marker='o', s=40, c='b', label = 'PpR3')
ax3.set_ylim([0,40])
ax3.set_xlim([0,20])
ax3.grid()
plt.legend(prop={'size': 14})


ax4=fig.add_subplot(344)
ax4.scatter(Del_R2_PpR3, ratio_phi_kex_PpR3, marker='o', s=40, c='b', label = 'PpR3')
ax4.set_ylim([0,40])
ax4.set_xlim([-2,20])
ax4.grid()
plt.legend(prop={'size': 14})

###############

ax5=fig.add_subplot(345)
ax5.scatter(kex_opt_PaR3, ratio_phi_kex_PaR3, marker='o', s=40, c='g', label = 'PaR3')
ax5.set_ylabel('$\Phi/K_{ex}$')
ax5.set_ylim([0,40])
ax5.set_xlim([-2,30])
ax5.grid()
plt.legend(prop={'size': 14})


ax6=fig.add_subplot(346)
ax6.scatter(Phiex_opt_PaR3, ratio_phi_kex_PaR3, marker='o', s=40, c='g', label = 'PaR3')
ax6.set_ylim([0,40])
ax6.set_xlim([-10,300])
ax6.grid()
plt.legend(prop={'size': 14})

ax7=fig.add_subplot(347)
ax7.scatter(R2o_opt_PaR3, ratio_phi_kex_PaR3, marker='o', s=40, c='g', label = 'PaR3')
ax7.set_ylim([0,40])
ax7.set_xlim([0,20])
ax7.grid()
plt.legend(prop={'size': 14})


ax8=fig.add_subplot(348)
ax8.scatter(Del_R2_PaR3, ratio_phi_kex_PaR3, marker='o', s=40, c='g', label = 'PaR3')
ax8.set_ylim([0,40])
ax8.set_xlim([-2,20])
ax8.grid()
plt.legend(prop={'size': 14})

###############

ax9=fig.add_subplot(3,4,9)
ax9.scatter(kex_opt_MfR3, ratio_phi_kex_MfR3, marker='o', s=40, c='r', label = 'MfR3')
ax9.set_xlabel('$K_{ex}$')
ax9.set_ylabel('$\Phi/K_{ex}$')
ax9.set_ylim([0,40])
ax9.set_xlim([-2,30])
ax9.grid()
plt.legend(prop={'size': 14})


ax10=fig.add_subplot(3,4,10)
ax10.scatter(Phiex_opt_MfR3, ratio_phi_kex_MfR3, marker='o', s=40, c='r', label = 'MfR3')
ax10.set_xlabel('$\Phi_{ex}$')
ax10.set_ylim([0,40])
ax10.set_xlim([-10,300])
ax10.grid()
plt.legend(prop={'size': 14})

ax11=fig.add_subplot(3,4,11)
ax11.scatter(R2o_opt_MfR3, ratio_phi_kex_MfR3, marker='o', s=40, c='r', label = 'MfR3')
ax11.set_xlabel('$R_{20}$')
ax11.set_ylim([0,40])
ax11.set_xlim([0,20])
ax11.grid()
plt.legend(prop={'size': 14})


ax12=fig.add_subplot(3,4,12)
ax12.scatter(Del_R2_MfR3, ratio_phi_kex_MfR3, marker='o', s=40, c='r', label = 'MfR3')
ax12.set_xlabel('$\Delta_{R2}$')
ax12.set_ylim([0,40])
ax12.set_xlim([-2,20])
ax12.grid()
plt.legend(prop={'size': 14})

plt.show()
############################


