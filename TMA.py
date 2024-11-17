# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 20:50:45 2024

@author: Larrot
"""
import numpy as np

def Tridiagonal_solver(a, b, c, d):
    """
     Algorithm to solve tridiagonal systems of equations

     Parameters
     ----------
     a : array
         Lower diagonal.
     b : array
         Middle diagonal.
     c : array
         Upper diagonal.
     d : array
         Right side of equations.

     Returns solution array x
     -------
    """
    size = len(d)
    x = np.zeros(size)
    # Forward steps
    for i in range(size-1):
        w = a[i]/b[i]
        b[i+1] -= w*c[i]
        d[i+1] -= w*d[i]
    # Backward steps
    x[size-1] = d[size-1]/b[size-1]
    for i in range(size-2, -1, -1):
        x[i] = (d[i]-c[i]*x[i+1])/b[i]

    return x