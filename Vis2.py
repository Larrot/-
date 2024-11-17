# # -*- coding: utf-8 -*-
# """
# Created on Fri Nov  1 18:59:14 2024

# @author: Larrot
# """

# # -*- coding: utf-8 -*-
# """
# Created on Fri Oct  4 21:04:25 2024

# @author: Larrot
# """




import CONSTS as C
import numpy as np
from Ionization_model import alpha
from Bisection import my_bisection
from scipy import optimize

# from scipy.optimize import fsolve





C_0 = C.K*C.N_A/C.MU*(27/(64*C.sigma*C.G*C.M_sun))**(1/4)

C_1 = (64*C.sigma*C.MU)/(27*C.K*C.N_A*(C.G*C.M_sun)**(1/2))




def T_diff(u, x, t):
    # nu_0 = 100
    # A = C.K*C.N_A/(C.MU*(C.G*C.M_sun)**(1/2))
    # def f(nu):
    # def f():
    #     #Example values:
    #     # Radial distance
    #     # x = 10
    #     # Surface density
    #     # u = 100
    #     # Time
    #     # t = 0.1
    #     return nu - alpha(u, x, t, nu)*T(u, x, t, nu)*(x*C.L)**(3/2)*A
    # f = lambda tu: tu - C_1*((x*C.L)**(3/2))/(C.Kappa*u**2*alpha1(u, x, t, nu))
    # f = lambda tu: tu - C_1*((x*C.L)**(3/2))/(C.Kappa*u**2*alpha(u, x, t, tu))*tu**4
    
    # f = lambda tu: tu - C_1*((x*C.L)**(3/2))/(C.Kappa*u**2*alpha(u, x, t, tu))*tu**4
    
    # f = lambda tu: 1 - C_1*((x*C.L)**(3/2))/(C.Kappa*u**2*0.01)*tu**3
    
    # f = lambda tu: tu - C_1*((x*C.L)**(3/2))/(C.Kappa*u**2*alpha(u, x, t, tu))*tu**4
    
    
    # f = lambda tu: tu - (C_1*((x*C.L)**(3/2))/(C.Kappa*u**2*0.01))**(-1/3)*tu**4
    f = lambda tu: tu**3 - 1/C_1*C.Kappa*u**2*alpha(u, x, t, tu)[0]*(x*C.L)**(-3/2)
    # def f(tu): return tu - 280*x**(-1/2)
    # f = lambda nu: nu - alpha(u, x, t, nu)*T(u, x, t, nu)*(x)**(3/2)
    # return alpha(u, x, t, nu)*T(u, x, t, nu)*(x*C.L)**(3/2)*A
    # solution = fsolve(f, nu_0)[0]
    # solution = my_bisection(f, 10**-8, 10**8, 0.001)
    solution = optimize.bisect(f, 10**-8, 10**8)
    return solution


def Diff(u, x, t):
    # D_0 = 5.11*10**(12)
    # D_0 = 1.83*10**(14)
    # T = T_diff(u, x, t)
    T = T_diff(u, x, t) + 280*x**(-1/2)
    # T = 280*x**(-1/2)
    # T = 17.45295*u**(2/3)*x**(-1/2)
    # T = (C_1*((x*C.L)**(3/2))/(C.Kappa*u**2*0.01))**(-1/3)
    # return (alpha(u, x, t, T)*C.N_A*C.K/(C.MU*(C.G*C.M_sun)**(1/2))*T*(x*C.L)**(3/2))
    # return (alpha(u, x, t, T)[0]*C.N_A*C.K/(C.MU*(C.G*C.M_sun)**(1/2))*T*(x*C.L)**(3/2))/D_0
    # return (alpha(u, x, t, T)[0]*C.N_A*C.K/(C.MU*(C.G*C.M_sun)**(1/2))*T*(x*C.L)**(3/2))/D_0
    return (alpha(u, x, t, T)[0]/0.0001)**(4/3)*(u/100)**(2/3)*x
    # return (0.01*T*x**(3/2))
    # return alpha(u, x, t, T)/0.0001*T*x**(3/2)/280
    # return 0.01/0.0001*T*x**(3/2)/280
    # return alpha(u, x, t, T)/0.0001*T*x**(3/2)/280
    # return alpha(u, x, t, T)/0.0001*x**(2/2)
# def Diff(u, x, t):

#     solution = (C_0*0.01)**(4/3)*(C.Kappa*u*(x*C.L))**(1/3)
#     return solution


# # Zero approximation
# # x0 = 100
# # Solve the equation
# # solution = fsolve(f, x0)
# # # print("Solution:", solution)
# DD_00 = C.N_A*C.K/(C.MU*(C.G*C.M_sun)**(1/2))*C.L**(3/2)
# print(f'{DD_00:.2e}')
# u = 100
# x = 100
# t = 0
# # # nu = 1000
# print(T_diff(u, x, t))

# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 19:37:39 2024

@author: Larrot
"""


# print(Diff(10**-5, 100, 0))






# TTT = (C_1*((x*C.L)**(3/2))/(C.Kappa*u**2*0.01))**(-1/3)
# print(TTT)
# print(T_diff(u, x, t))
# # print(Diff(u, x, t, TTT))
# # print(Diff(u, x, t, T_diff(u, x, t)))
# # print(Ionization_Rate(u, x, t, T_diff(u, x, t)))
# # nu = Diff(u, x, t)
# # nu = 1000
# # print(C_0*alpha(u, x, t, nu)*(C.Kappa*nu*u*(x*C.L))**(1/4))
# # print((C_0*0.01)**(4/3)*(C.Kappa*u*(x*C.L))**(1/3))
# # nu = (C_0*0.01)**(4/3)*(C.Kappa*u*(x*C.L))**(1/3)

# # print(T(u, x, t, nu))
# # print(Ionization_Rate(u, x, t, nu))


