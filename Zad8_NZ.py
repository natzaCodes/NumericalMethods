"""
Created on Fri May 13 14:40:24 2022
@author: Natalia
Kwadratury Gaussa-Legendre'a
"""
import numpy as np
import scipy.special



def GL(f,a,b):
    n=5
    x,w=scipy.special.roots_legendre(n) #zwraca węzły xi oraz wi - to co było w zadaniu domowym
    I=0
    for i in range(n):
        I+=w[i]*f((b-a)*x[i]/2.+(a+b)/2.)
    return I*(b-a)/2.
    
def GL_rec(f,a,b,eps,I0):
    L=GL(f,a,(a+b)/2.)
    P=GL(f,(a+b)/2.,b)
    
    if np.abs(L+P-I0)/np.abs(I0+1e-15) < eps:
        return L+P
    else:
        I=GL_rec(f,a,(a+b)/2.,eps,L)+GL_rec(f,(a+b)/2.,b,eps,P)
        return I
    
    

def f(x):
    return np.cos(x**2)

I0=GL(f,-1,1)
I=GL_rec(f,-1,1,1e-5,I0)
print(I0)
print(I)