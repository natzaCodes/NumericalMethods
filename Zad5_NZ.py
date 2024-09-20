"""
Created on Mon Apr 18 16:23:59 2022
@author: Natalia
Zadanie zaliczeniowe 5 - wyznaczanie miejsc zerowych
"""

import numpy as np
import sys

'''
 można uniknąć wpadania w nieskończone pętle np. poprzez ustalenie maksymalnej liczby iteracji, która będzie duża, 
 ale nie za duża :)
'''

def bisekcja(f, a, b, eps):
    x=12.5 #potrzebuje czegos bardzo randomowego
    y=12.5 #potrzebuje czegos bardzo randomowego
    s=0
    s_max = 1000
    c2 = (a+b)/2
    xR = 42282.703767501494
    if f(a)*f(c2)>0 and f(c2)*f(b)>0:  
        print("Przedział jest niepoprawny. Wartoci funkcji w granicach przedziału muszą mieć przeciwne znaki - met. bisekcji.")
        sys.exit()
    
    while abs(b-a)/2. > eps and s<s_max:  
        s+=1
        c=(a+b)/2    
        if f(c)==0:
            x,y = c,f(c)
            break
        elif f(a)*f(c)<0:
            b=c  
            x,y = c,f(c)
        else:
            a=c
            x,y = c,f(c)
     
    eps2=abs(xR-x)
    return x,y,s,eps2


def falsi(f, a, b, eps):
    x=12.5 #potrzebuje czegos bardzo randomowego
    y=12.5 #potrzebuje czegos bardzo randomowego
    s=0
    s_max = 100
    xR = 42282.703767501494
    
    if f(a)*f(b)>0 :  
        print("Przedział jest niepoprawny. Wartoci funkcji w granicach przedziału muszą mieć przeciwne znaki - reguła falsi.")
        sys.exit() 
    
    while abs(a-b) > eps:
        if s==s_max:
            break
        s+=1
        #y0=a0+a1*x
        a1 = (f(b)-f(a))/(b-a)
        a0 = f(a)-a1*a
        #pierwsze przybliżenie pierwiastka
        x1=(0-a0)/a1 
        
        if f(x1)==0:
            x,y = x1,f(x1)
            break
        elif f(a)*f(x1)<0:
            b=x1  
            x,y = x1,f(x1)
        else:
            a=x1
            x,y = x1,f(x1)
    
        eps2=abs(xR-x)
    return x,y,s,eps2

def sieczne(f, a, b, eps):
          
    x=12.5 #potrzebuje czegos bardzo randomowego
    y=12.5 #potrzebuje czegos bardzo randomowego   
    #print("xR=%.15f" % xR)
    x0=a
    x1=b
    s=0 
    s_max = 500
    while abs(x0-x1) > eps: 
        if s==s_max:
                break
        s+=1
        x = x1-f(x1)*(x1-x0)/(f(x1)-f(x0)+1e-20) #ryzyko dzielenia przez zero
        if f(x1)-f(x0) ==0:
            print("nastąpiło dzielenie przez 0 w metodzie siecznych w kroku",s)
        y=f(x)
        x0 = x1
        x1 = x
        eps2 = abs(x0-x1)
        
        #print("x=", round(x,8), "epsilon",s," =", round(eps2, 8))

    return x,y,s,eps2

def Newton(f, df, x0, eps):
    
    x=x0 
    y=f(x0) 
    s=0
    s_max = 1000
    x1 = x0+1
    while abs(x-x1) > eps:
        if s==s_max:
                break
        s+=1
        x1=x
        if f(x0) == 0:
            x,y = x0, f(x0)
        else:
            x = x-f(x)/(df(x)+1e-20)  #ryzyko dzielenia przez zero
            if df(x) ==0:
                print("nastąpiło dzielenie przez 0 w metodzie Newtona w kroku",s)
            y=f(x)
        eps2 = abs(x1-x)
        
    
    return x,y,s,eps2

#zdefiniowana funkcja:
def f(x):
    mv=15.204
    Mv= -3.093
    alphav= 3.93036e-7
    return x-np.power(10,(mv-Mv+5-alphav*x)/5)

def df(x):
    mv=15.204
    Mv= -3.093
    alphav= 3.93036e-7
    return 1+np.power(10,(mv-Mv+5-alphav*x)/5)*np.log(10)*alphav/5


#do testów
a=1e4
b=1e5
eps=1e-10
x0=1e4

print("x, y, l. kroków, błąd:",bisekcja(f, a, b, eps))
print("x, y, l. kroków, błąd:",falsi(f, a, b, eps))
print("x, y, l. kroków, błąd:",sieczne(f, a, b, eps))
print("x, y, l. kroków, błąd:",Newton(f, df, x0, eps))


'''
Metoda bisekcji wymaga najwiekszej liczby kroków do obliczenia miejsca zerowego oraz posiada
najmniejszą dokładnosć wraz z metodą falsi. Metody siecznnych i stycznych sa najbardziej dokładne.
'''










