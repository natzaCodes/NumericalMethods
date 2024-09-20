"""
Created on Fri May 27 2022
@author: Natalia Zalewska
ordinary differential equations solver
"""
import numpy as np


def Rk1(f,h,t,x0,y0,z0):
    n = len(t)
    for i in range(n-1):  
        y= y0+h*f(x0,y0)
        y0=y
        x0+=h
    return y
    
def Rk2(f,h,t,x0,y0,z0):
    n = len(t)
    
    for i in range(n-1):
        k1= f(x0,y0)
        y=y0+h*f((x0+h/2), (y0+h/2*k1))
        y0=y
        x0+=h
    return y

def Rk3(f,h,t,x0,y0,z0):
    n = len(t)
    
    for i in range(n-1):
        k1= f(x0,y0)
        k2 = f((x0+h/2), (y0+h/2*k1))
        k3 = f((x0+h*3/4), (y0+h*3/4*k2))
        k=(2*k1+3*k2+4*k3)/9  
        y=y0+h*k
        y0=y
        x0+=h
    return y

def Rk4(f,h,t,x0,y0,z0):
    n = len(t)
    
    for i in range(n-1):
        k1 = f(x0, y0)
        k2 = f((x0+h/2), (y0+k1*h/2))
        k3 = f((x0+h/2), (y0+k2*h/2))
        k4 = f((x0+h), (y0+h*k3))
        k = (k1+2*k2+2*k3+k4)/6
        y = y0 + h*k
        y0 = y
        x0 = x0+h
    return y
        
def f(x,y,z):
    return
def z(x,y):
    return    

h=1e-4
ksi=np.arange(0,10,h)































