"""
Created on Fri Apr  8  2022
Neville's algorithm, drawing a given and interpolated polynomial
@author: Natalia Zalewska
"""

import numpy as np
import matplotlib.pyplot as plt

def Neville (wn,fn,x):
    p = fn.copy()
    n = len(p)

    for j in range(1,n):
        for i in range(n-j):
            p[i] = ((wn[i+j] - x)* p[i] +(x- wn[i])*p[i+1]) / (wn[i+j] - wn[i])
    return p[0]

#print(Neville(w,f,0))

def Lagrange(wn, fn):
    sortix = np.argsort(np.abs(wn))[::-1]
    wn, p = wn[sortix], fn[sortix]
    
    #sortowanie tablicy wn i fn wzgledem modułu wn a potem
    
    n=len(wn)
    a = np.zeros(n)
    
    a[0]=Neville(wn,p, x=0)
    
    
    for j in range(1,n):
        for i in range(n-j):
            p[i] = (p[i] - a[j-1])/ wn[i]
        a[j] = Neville(wn[:n-j],p[:n-j],x=0)
            
    return a

def polyncalc(x,N,a):  
    #schemat Hornera
    n=N
    P=0
    while n>=0:
        P = P*x + a[n]
        n-=1
    return P

path = "dane.txt"
w, f = np.loadtxt(path, unpack=True)
N =len(w)-1

''' podpunkt a '''

n=0
a10 = []
while n<=N:
    print("Podaj", n, "współczynnik dla wielomianu", N, "stopnia:")
    a = float(input())
    a10.append(a)
    n+=1

#policzenie wartosci dla x z pliku zadanego wielomianu
Ptab1 = np.ones(N+1)
n=0
for x in w:
    Ptab1[n] = polyncalc(x,N,a10)  
    n+=1
    
#liczenie nowych współczynników wielomianu
a1 = Lagrange(w, Ptab1)

xmin1 = min(w)+10*min(w)
xmax1 = max(w)+10*max(w)
kroki = int((xmax1 - xmin1)*100)
X1=np.linspace(xmin1, xmax1, kroki)
Ptab11 = [] 
for x in X1:
    P=polyncalc(x,N,a1)
    Ptab11.append(P)


print("\n Stare współczynniki to:", a10)
print("\n Nowe współczynniki to:", a1)



''' podpunkt b '''

xmin = -1
xmax = 1
step = (xmax-xmin)/(N)
X=np.arange(xmin, xmax+0.2, step)

#wspolczynniki wielomianu dla sin dla dowolnego stopnia wielomianu
a20 = np.zeros(N+1)
a20[0]=0
a20[1]=1
s=0
for n in range(2,N+1):
    if n%2==0:
        a20[n]=0
    else: 
        s+=1
        a20[n] = (-1)*a20[n-2]/(2*s*(2*s+1))
        
    
#policzenie wartosci
Ptab2 = np.ones(N+1)
n=0
for x in X:
    Ptab2[n] = polyncalc(x,N,a20)  
    n+=1

    
#liczenie nowych współczynników wielomianu
a2 = Lagrange(X, Ptab2)
kroki = int((xmax - xmin)*100)
X2=np.linspace(xmin-6, xmax+6, kroki)
Ptab22 = [] 
for x in X2:
    P=polyncalc(x,N,a2)    
    Ptab22.append(P)

print("\n Stare współczynniki dla sin to:", a20)
print("\n Nowe współczynniki dla sin to:", a2)


plt.figure()
plt.scatter(w, Ptab1, s=10, label='Wprowadzone współczynniki')
plt.plot(X1, Ptab11, color='maroon', label='Po interpolacji')
plt.title("Wielomian wprowadzony ręcznie",fontsize=10, fontweight=0) 
plt.xlabel('x')
plt.ylabel('Wartosc')  
plt.legend()
plt.savefig("Wykres_NZ21.png",bbox_inches='tight', pad_inches=0) 
plt.draw()


plt.figure()
plt.scatter(X, Ptab2,s=20, label='Z Taylora')
plt.plot(X2, Ptab22, color='maroon', label='Po interpolacji')
plt.title("Funkcja sinus jako wielomian",fontsize=10, fontweight=0) 
plt.xlabel('x')
plt.ylabel('Wartosc')  
plt.legend()
plt.savefig("Wykres_NZ2.png",bbox_inches='tight', pad_inches=0) 
plt.show()






'''
W punkcie a współczynniki wychodzą takie same (lub bardzo zbliżone), w podpunkcie b wykres jest ładny (tam, gdzie przechodzi przez punkty), 
współczynniki w sumie są bardzo zbliżone
'''













