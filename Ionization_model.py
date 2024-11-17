import CONSTS as C

import numpy as np


def sigmoid(x, x0, k):

    return 1 / (1 + np.exp(-k * (x - x0)))

def S_ALPHA(x):
    y2 = 0.01
    y1 = 0.0001
    transition_point = -12 
    # transition_point = 10**(-1) 
    # steepness =   
    steepness = 20
    
    return y1 + (y2 - y1) * sigmoid(x, transition_point, steepness)


def alpha(u, x, t, T):
    Ksi = C.Ksi_0*np.exp(-u/C.R_CR)

    n = C.C1*u*T**(-1/2)*x**(-3/2)

    if 0 <= T <= 10**3:
        alpha_r = 2.07*10**(-11)*T**(-1/2)*3
    elif 10**3 <= T <= 10**4:
        alpha_r = 2.07*10**(-11)*T**(-1/2)*1.5
    else:
        alpha_r = 0

    if T <= 150:
        alpha_g = C.alpha_g0

    elif 150 <= T <= 400:
        alpha_g = -1/250*(C.alpha_g0-C.alpha_gm)*T + \
            8/5*C.alpha_g0 - 3/5*C.alpha_gm

    elif 400 <= T <= 1500:
        alpha_g = C.alpha_gm

    elif 1500 <= T <= 2000:
        alpha_g = -1/500*C.alpha_gm*T+4*C.alpha_gm
    else:
        alpha_g = 0

    # betta
    if alpha_r != 0:
        betta = (alpha_g*n+Ksi) / (2*alpha_r*n)
        betta1 = (alpha_g*n+C.Ksi_r) / (2*alpha_r*n)
        gama = Ksi/(alpha_r*n)
        gama1 = C.Ksi_r/(alpha_r*n)
    else:
        betta = 0
        betta1 = 0
        gama = 0
        gama1 = 0

    Xi = -betta + (betta**2+gama)**(1/2)
    Xi_r = -betta1 + (betta1**2+gama1)**(1/2)
    Xt = 1.8*10**(-11)*(T/1000)**(3/4)*(C.Nu_K/10**(-7))**(0.5) \
        * (n/10**(13))**(-0.5) \
        * np.exp(-25000/T)/(1.15*10**(-11))

    X = Xi + Xi_r + Xt
    
    if X == 0 or X >= 1:
        X = 1

    if X < 10**(-12):
        ALPHA = 10**(-4)
        # ALPHA = 10**(-2)
        # ALPHA = 0.1*10**(-2)
    else:
        ALPHA = 0.01
        # ALPHA = 10**-4
    return [ALPHA, X]
    # return [S_ALPHA(np.log10(X)), X]
