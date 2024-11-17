# # -*- coding: utf-8 -*-
# """
# Created on Mon Nov  4 22:17:37 2024

# @author: Larrot
# """
import CONSTS as C
import numpy as np
from Vis2 import T_diff, Diff
from Ionization_model import alpha

t1 = 0

t2 = 1

t3 = 2

t4 = 3


T = np.array([C.tau*i for i in range(C.M)])

TAU = 1.39


DATA = np.loadtxt(f"DATA_{C.t_end}_{C.M}_{C.Dead}.txt")
def u_init(x):
    if x <= 100:
        u_init = 100*x**(-3/8)
        # u_init = 10**(-3/8)*100*x**(-3/8)
        # u_init = 100*x**(-3/8)
        # u_init = 100
    else:
        u_init = 10**(-5)
        # u_init = 0
    return u_init



x_s = np.logspace(C.R_in, C.R_out, C.N+1)
# x = np.logspace(R_in, R_out, N+1)
x = [(x_s[i+1]+x_s[i])/2 for i in range(C.N)]

# DATA = np.zeros((C.M, C.N))

# for t in range(C.M):
#     for i in range(C.N):
#         DATA[t][i] = u_init(x_s[i])

TT = np.zeros((C.M, C.N))

for t in range(C.M):
    for i in range(C.N):
        TT[t][i] = T_diff(DATA[t][i], x_s[i], t)
        
ALPHA = np.zeros((C.M, C.N))
for t in range(C.M):
    for i in range(C.N):
        ALPHA[t][i] = alpha(DATA[t][i], x_s[i], 0, TT[t][i])[0]
    
XX = np.zeros((C.M, C.N))
for t in range(C.M):
    for i in range(C.N):
        XX[t][i] = alpha(DATA[t][i], x_s[i], 0, TT[t][i])[1]
    
DIFF = np.zeros((C.M, C.N))
for t in range(C.M):
    for i in range(C.N):
        DIFF[t][i] = Diff(DATA[t][i], x_s[i], 0)

from matplotlib import pyplot as plt
import colorcet as cc

clrs = []
for i in range(4):
    """ 
    выбираем палитру fire
    
    всего в палитре 256 цветов, 
    поэтому если хотим сделать список из 3х цветов из этой палитры, 
    то 256/3 ~ 64 - с таким шагом надо выбрать из палитры цвета
    """
    clrs.append(cc.fire[::][i*64])

colors = [clrs[3], clrs[2], clrs[1], clrs[0]]

    
fig = plt.figure(figsize=(8, 6))
# fig = plt.figure(figsize=(10, 8))
# fig = plt.figure(figsize=(6, 4))
# fig.subplots_adjust(hspace=0.6, wspace=0.4)

# """График поверхностной плотности"""
plt.rcParams.update({'font.size': 17})

# ax_11 = fig.add_subplot(1, 1, 1)
# ax_11.plot(x[1:], TT[t1][1:], '-', color=colors[1],
#           label=fr'$t=${TAU*T[t1]:.2f} млн лет')
# ax_11.plot(x[1:], TT[t2][1:], '-.', color=colors[2],
#           label=fr'$t=${TAU*T[t2]:.2f} млн лет')
# ax_11.plot(x[1:], TT[t3][1:], '--', color=colors[3],
#           label=fr'$t=${TAU*T[t3]:.2f} млн лет')
# ax_11.plot(x[1:], TT[t4][1:], '-', color=colors[3],
#           label=fr'$t=${TAU*T[t4]:.2f} млн лет')


# ax_11.plot(x[1:], B_z4[t1][1:], '--', color='k',
#           label=f'$B_z$ для min')

# ax_11.set_ylim(10**(-6), 10**(6))
# ax_11.set_xlim(0.01, 1000)


# ax_11.set_xlabel('$r,~{а.е.}$')
# ax_11.set_ylabel(r'$T,~{К}$')

# ax_11.set_yscale('log')
# ax_11.set_xscale('log')

# ax_11.set_title(fr'Температура')
# ax_11.legend(loc='best')

ax_11 = fig.add_subplot(1, 1, 1)
ax_11.plot(x[1:], XX[t1][1:], '-', color=colors[1],
          label=fr'$t=${TAU*T[t1]:.2f} млн лет')
ax_11.plot(x[1:], XX[t2][1:], '-.', color=colors[2],
          label=fr'$t=${TAU*T[t2]:.2f} млн лет')
ax_11.plot(x[1:], XX[t3][1:], '--', color=colors[3],
          label=fr'$t=${TAU*T[t3]:.2f} млн лет')
ax_11.plot(x[1:], XX[t4][1:], '-', color=colors[3],
          label=fr'$t=${TAU*T[t4]:.2f} млн лет')


# # # ax_11.plot(x[1:], B_z4[t1][1:], '--', color='k',
# # #           label=f'$B_z$ для min')

# # # ax_11.set_ylim(10**(-6), 10**(6))
ax_11.set_xlim(0.01, 1000)


ax_11.set_xlabel('$r,~{а.е.}$')
ax_11.set_ylabel(r'$x$')

ax_11.set_yscale('log')
ax_11.set_xscale('log')

ax_11.set_title(fr'Степень ионизации')
ax_11.legend(loc='best')

# ax_11 = fig.add_subplot(1, 1, 1)
# ax_11.plot(x[1:], ALPHA[0][1:], '-', color=colors[1],
#            label=fr'0')
# ax_11.plot(x[1:], ALPHA[5][1:], '-.', color=colors[2],
#           label=fr'5')
# ax_11.plot(x[1:], ALPHA[10][1:], '--', color=colors[3],
#           label=fr'10')


# ax_11.plot(x[1:], B_z4[t1][1:], '--', color='k',
#           label=f'$B_z$ для min')

# ax_11.set_ylim(10**(-6), 10**(6))
# ax_11.set_xlim(0.01, 1000)


# ax_11.set_xlabel('$r,~{а.е.}$')
# ax_11.set_ylabel(r'$\alpha$')

# ax_11.set_yscale('log')
# ax_11.set_xscale('log')

# ax_11.set_title(fr'Коэфф. альфа dead zone')
# ax_11.legend(loc='best')

fig.savefig(f"X.png", orientation='landscape', dpi=300)

