"""
Created on Fri Apr  1  2022
interpolation polynomial - points with the polynomial on the graph
@author: Natalia Zalewska
"""

import numpy as np
import sys
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

#print(Neville(w,f,0))

def polyncalc(w,f,x,N):  
    #schemat Hornera
    n=N
    P=0
    while n>=0:
        P = P*x + Lagrange(w,f)[n]
        n-=1
    return P



''' Main code '''

path = "dane.txt"
w, f = np.loadtxt(path, unpack=True)
N =len(w)-1

if len(w)>6:
    print("Jest ok, lecimy dalej")
else:
    print("Werszy ma być więcej niż 6, podaj inny plik" )
    sys.exit()



print("Współczynniki wielomianu interpolacyjnego to:")
i=0
while i < len(w):
    print(i,".",round(Lagrange(w,f)[i], 4))
    i+=1

i=1
x0t=[]
Pt=[]
while i>0:    
    print("\n Podaj x dla którego chcesz poznać wartosć wielomianu interpolującego:")
    x0 = np.float32(input())
    P = polyncalc(w,f,x0,N)
    x0t.append(x0)
    Pt.append(P)
    print("\n Wartosc wielomianu dla x=:",x0, "to:",P)  
    print("\n Czy chcesz sprawdzić kolejny punkt? \n 1=TAK, 0=NIE")
    i = np.int(input())    
 
xmin = min(min(w),min(x0t))
xmax = max(max(w),max(x0t))
kroki = int((xmax - xmin)*100)
X=np.linspace(xmin, xmax, kroki)
Ptab = [] 

for x in X:
    P = polyncalc(w,f,x,N)  
    Ptab.append(P)


# przyda sie do wykresu
napis = "Wspolczynniki to: \n"
n=0
while n<=N:
    napis = napis + str(round(Lagrange(w,f)[n],2)) + ', '
    n+=1
    
    
    
plt.plot(X, Ptab, label='Dopasowany wielomian')
plt.scatter(w, f, s=20, color='maroon', label='Wczytane punkty')
plt.scatter(x0t, Pt, s=20, color='darkgreen', label='Punkty podane przez użytkownika')
plt.title("Wielomian interpolacyjny stopnia {0}. \n {1} \n".format(N, napis),fontsize=10, fontweight=0) 
plt.xlabel('x')
plt.ylabel('Wartosc')  
plt.legend()
plt.savefig("Wykres1_NZ.png",bbox_inches='tight', pad_inches=0) 
















