import numpy as np
import pylab
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import FixedLocator
import string
pylab.rcParams['figure.figsize'] = (10.0, 5.0)

# Parameterization
K_max = 10000
r_max = 0.45
d = 0.01
g = 0.8
sigma  = 1
k = 2
b = 10
aRS = 0.9
aSR = 0.1
# Initialization
u = np.linspace(0,1,1001)

# Ecological equilibrium points for (m,u) values in [0,1]*[0,1]
def f(Y, K_max, r_max, d, g, b, k):
    v1, m1 = Y
    x_R_temp = max(0,(K_max/(r_max*(1-aSR*aRS)))*(aSR*m1*np.exp(g*v1)/(k+b*v1)\
             - m1/k + aSR*d*np.exp(g*v1) + r_max -d - aSR*r_max))
    x_S_temp = max(0,(K_max/(r_max*(1-aSR*aRS)))*(-m1*np.exp(g*v1)/(k+b*v1) \
            + aRS*m1/k - d*np.exp(g*v1) + r_max - aRS*r_max + aRS*d))
    return [x_R_temp + x_S_temp]

# Meshgrid to determine ecological equilibrium values
v = np.linspace(0, 1, 300)
m = np.linspace(0, 1.005, 300)
Y1, Y2 = np.meshgrid(m, v)
NI, NJ = Y1.shape
x = np.zeros(np.shape(Y1))
for i in range(NI):
    for j in range(NJ):
        v_m = Y2[i, j]
        m_m = Y1[i, j]
        yprime = f([v_m, m_m], K_max, r_max, d, g, b, k)
        x[i,j] = yprime[0]
# Progression, Stabilization, and Extinction regions
z = x
levels = FixedLocator([z.min(), 0, 7000, z.max()], nbins=4).tick_values(z.min(), z.max())
colorsList = ['green', '#FAFA52', 'tomato']
cmap = matplotlib.colors.ListedColormap(colorsList)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
im = ax1.pcolormesh(Y1, Y2, z, cmap=cmap, norm=norm, shading='nearest')
im = ax2.pcolormesh(Y1, Y2, z, cmap=cmap, norm=norm, shading='nearest')
#Followers's best response curve
mm1 = np.linspace(0,1,1001)
u_Stack_temp = -k/b - mm1/(2*b*d) + np.sqrt((mm1 * g)**2 + 4*mm1*g*b*d)/(2*b*d*g)
# Illustrating Follower's best response curve
ax1.plot(mm1,u_Stack_temp,color='b', linewidth=3)
ax2.plot(mm1,u_Stack_temp,color='b',label="Followers' Best Response", linewidth=3)


#####################################################################
######################### QoL with uR, c2 > 0 #######################

# Quality of life coefficients
c1 = 0.54
c2 = 0.21
c3 = 0.25
# Leader's best response curve
kk = (1/(1-aSR*aRS))**2
r_u = r_max*np.exp(-g*u)
m = (c1*kk*((aSR-1)/((k+b*u)*r_u)+(aRS-1)/(k*r_max))*\
    ((1-aRS)*(1-d/r_max)+(1-aSR)*(1-d/r_u)))/\
        (-c1*kk*((aSR-1)/((k+b*u)*r_u)+(aRS-1)/(k*r_max))*\
            ((aRS-1)/(k*r_max)+(aSR-1)/((k+b*u)*r_u))-c3)
# Figure 4A, Leader's best response curve illustration
ax1.plot(m,u,'b--', linewidth=3)
# MTD solution index, m=1
m_MTD_idx = 1000
# Stackelberg solution
m_st_idx = 660
u_st = u_Stack_temp[m_st_idx]
# x* value at Stackelberg solution
x_st = np.array(f([u_st,mm1[m_st_idx]], K_max, r_max, d, g, b, k))
# QoL at Stackelberg solution
Q_S = float(1 - c1*((x_st/K_max)**2) - c2*u_st**2 - c3*mm1[m_st_idx]**2)
# Illustration
ax1.scatter(mm1[m_st_idx],u_st,color='b', s=50)
ax1.text(mm1[m_st_idx]-0.182,u_st+0.06,'Stackelberg',backgroundcolor='1',size='small')
ax1.text(mm1[m_st_idx]-0.183,u_st+0.018,'Q='+str(round(Q_S,3)),backgroundcolor='1',size='x-small')
# Nash solution
m_Na_idx = 750
u_Na = u_Stack_temp[m_Na_idx]
# x* value at Nash solution
x_Na = np.array(f([u_Na,mm1[m_Na_idx]], K_max, r_max, d, g, b, k))
# QoL at Nash solution
Q_Na = float(1 - c1*((x_Na/K_max)**2) - c2*u_Na**2 - c3*mm1[m_Na_idx]**2)
# Illustration
ax1.scatter(mm1[m_Na_idx],u_Na,color='b', s=50)
ax1.text(mm1[m_Na_idx]+0.028,u_Na-0.05,'Nash',backgroundcolor='1',size='small')
ax1.text(mm1[m_Na_idx]+0.028,u_Na-0.1,'Q='+str(round(Q_Na,3)),backgroundcolor='1',size='x-small')
ax1.scatter(mm1[m_MTD_idx],u_Stack_temp[m_MTD_idx],color='b', s=50)
ax1.text(mm1[m_MTD_idx]-0.09,u_Stack_temp[m_MTD_idx]+0.03,'MTD',backgroundcolor='1',size='small')


####################################################################
################## QoL without uR, c2 = 0 ##########################

# Quality of life coefficients
c1 = 0.68
c3 = 0.32
# Leader's best response curve
kk = (1/(1-aSR*aRS))**2
r_u = r_max*np.exp(-g*u)
m = (c1*kk*((aSR-1)/((k+b*u)*r_u)+(aRS-1)/(k*r_max))*\
    ((1-aRS)*(1-d/r_max)+(1-aSR)*(1-d/r_u)))/\
        (-c1*kk*((aSR-1)/((k+b*u)*r_u)+(aRS-1)/(k*r_max))*\
            ((aRS-1)/(k*r_max)+(aSR-1)/((k+b*u)*r_u))-c3)
# Figure 4B, Leader's best response curve illustration
ax2.plot(m,u,'b--',label="Leader's Best response", linewidth=3)
# Nash and Stackelberg solutions coincide
m_Nb_idx = 743
u_Nb = u_Stack_temp[m_Nb_idx]
# x* at Nash solution
x_Nb = np.array(f([u_Nb,mm1[m_Nb_idx]], K_max, r_max, d, g, b, k))
# QoL at Nash solution
Q_N = float(1 - c1*((x_Nb/K_max)**2) - c3*mm1[m_Nb_idx]**2)
# Illustration
ax2.scatter(mm1[m_Nb_idx],u_Nb,color='b', s=50)
ax2.text(mm1[m_Nb_idx]+0.028,u_Nb-0.05,'Nash',backgroundcolor='1',size='small')
ax2.text(mm1[m_Nb_idx]+0.028,u_Nb-0.1,'Stackelberg',backgroundcolor='1',size='small')
ax2.text(mm1[m_Nb_idx]+0.028,u_Nb-0.15,'Q='+str(round(Q_N,3)),backgroundcolor='1',size='x-small')
ax2.scatter(mm1[m_MTD_idx],u_Stack_temp[m_MTD_idx],color='b', s=50)
ax2.text(mm1[m_MTD_idx]-0.1,u_Stack_temp[m_MTD_idx]+0.03,'MTD',backgroundcolor='1')
ax2.set_xlim([0,1.01])
#######################################
################################
########################
ax1.set_xlim([0,1.01])
ax1.text(-0.1, 1.05, string.ascii_uppercase[0], transform=ax1.transAxes, 
            size=20, weight='normal')
ax2.text(-0.1, 1.05, string.ascii_uppercase[1], transform=ax2.transAxes, 
            size=20, weight='normal')
ax1.set_xlabel('Treatment dose (m)', size = 12)
ax1.set_ylabel('Resistance level (u)', size = 12)
ax1.set_ylim((0,1))
ax1.set_xlim((0,1))
ax2.set_xlabel('Treatment dose (m)', size = 12)
ax2.set_ylabel('Resistance level (u)', size = 12)
ax2.set_ylim((0,1))
ax2.set_xlim((0,1))
fig.legend()

plt.savefig("subplot_cancer.pdf", dpi=100)

plt.show()