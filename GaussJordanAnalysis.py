"""
Created on Fri Apr 29 2022
@author: Natalia Zalewska
Application of the Gauss Jordan method
"""
import numpy as np
import numpy.linalg
import sys
import matplotlib.pyplot as plt


def Gauss_Jordan(macierz_glowna, data,czy_liczyc_macierz_odwrtona):
    
    matrix=data  
    diagonalna = np.diag(np.ones(len(macierz_glowna[:,1])) )
    odwracanie = np.concatenate((macierz_glowna,diagonalna), axis=1)
    N = len(matrix)
    
 
    if czy_liczyc_macierz_odwrtona==False:
        for j in range(N):
            if matrix[j,j]==0:
                sys.exit()
            matrix[j,j:]/=matrix[j,j] 
            for k in range(0,N):
                if k!=j:
                    matrix[k,:]=matrix[k,:]-matrix[j,:]*matrix[k,j]
        return macierz_glowna, matrix[:,len(matrix[1,:])-1]
    
    if czy_liczyc_macierz_odwrtona==True:
        for j in range(N):
            if odwracanie[j,j]==0:
                sys.exit()
            odwracanie[j,j:]/=odwracanie[j,j]
            for k in range(0,N):
                if k!=j:
                    odwracanie[k,:]=odwracanie[k,:]-odwracanie[j,:]*odwracanie[k,j]
                    
        return macierz_glowna, matrix[:,len(matrix[1,:])-1],odwracanie[:,len(matrix):len(odwracanie[1,:])]
    
    

def MNK_polyn(m, N, x, y):
        
    Sx= np.zeros(2*m+1)
    Sxy=np.zeros((m+1,1))
    A=np.zeros((m+1,m+1))
    '''
    #to z jakeigos powodu wywala w kosmos
    x1=np.copy(x)
    Sx[0]=N
    Sx[1]=np.sum(x)
    Sxy[0,0]=np.sum(y)
    for k in range(1,m):
        for j in range(1,m+1):
            Sxy[j,0]=0
            Sx[k+j]=0
            for i in range(N-1):
                Sx[k+j]+= x[i]
                Sxy[j,0]+=x[i]*y[i]
            x*=x1 #liczy x do kolejnych potęg
    '''
    for j in range(m+1):
        for k in range(m+1): 
            Sx[k+j]= np.sum(np.power(x,k+j))
        Sxy[j,0]=np.sum(np.power(x,j)*y)
            
    for j in range(m+1):
        for k in range(m+1):
            A[j,k]=Sx[j+k]
    
    
    return Sx, Sxy, A

def polyn_write(N,atab,m):
    
    print("\n Twój wielomian stopnia",m," to:")
    n=0
    while n<=N:
        
        if atab[n]<0:
            print(" ", end="")
        elif atab[n] == 0 or n==0:
            print(end="") 
        else:
            print(' + ', end="")
            
        if atab[n] != 0 and atab[n] !=1 and atab[n] !=(-1) and n>=1:
            print(atab[n], end="")    
        elif atab[n] !=0 and n==0:
            print(atab[n], end="")
        elif atab[n] ==1 and n>=1:
            print(end="")
        elif atab[n] == (-1) and n>=1:
            print("-", end="")
            
       
        if n == 0 or atab[n] == 0:
            print(end="")
        elif n == 1:
            print('x', end="")
        else:
            print('x^', end="") 
            print(n, end="")
         
        n+=1
    
    return 

def Horner(m,X,a):
    
    Pm = []
    for x in X:
        n=m-1
        P=a[m]
        while n>=0:
            P = P*x + a[n]
            n-=1
        Pm.append(P)
    
    return Pm

def Statistics(X,Y,P):
    E= np.sum((Y-P)**2)
    RMSD=np.sqrt(E/len(X))
    return RMSD,E

'''Main code'''

# cm, bar, g/cm^3, K 
X, Y, ro, T = np.loadtxt("Zadanie_Zaliczeniowe_7.dat",unpack=True)
n=np.size(X)-1
m1=5
m2=10
X=X/1e10
Y=Y/1e11

#dla m=5
Sx1,Sxy1,A1=MNK_polyn(m1,n,X,Y)
matrix_gauss1 = np.concatenate((A1,Sxy1), axis=1)
A1g, a1 =Gauss_Jordan(A1, matrix_gauss1,False)
Sx1g= np.zeros(2*m1+1)
for w in range(m1+1) :
    for k in range(m1+1):
        Sx1g[w+k]=A1g[w,k]

#dla m=10
Sx2,Sxy2,A2=MNK_polyn(m2,n,X,Y)
matrix_gauss2 = np.concatenate((A2,Sxy2), axis=1)
A2g, a2 =Gauss_Jordan(A2, matrix_gauss2,False)
Sx2g= np.zeros(2*m2+1)
for w in range(m2+1) :
    for k in range(m2+1):
        Sx2g[w+k]=A2g[w,k]

#oblicznie wartoci g(x) dla X z tabeli na podstawie wyznaczonego wielomianu
P1 = Horner(m1,X,a1)
P2 = Horner(m2,X,a2)
RMSD1,E1 = Statistics(X,Y,P1)
RMSD2,E2 = Statistics(X,Y,P2)
 



print("\n Wartosci przed zastosowaniem algorytmu dla m=5 to: \n Sx = ", np.matrix.round(Sx1), "\n Sxy =", np.around(Sxy1))
print("\n Wartosci po zastosowaniu algorytmu to: \n Sx = ", np.matrix.round(Sx1g), "\n a =", np.around(a1,5))
print("\n Wartosci przed zastosowaniem algorytmu dla m=10 to: \n Sx = ", np.matrix.round(Sx2), "\n Sxy =", np.around(Sxy2))
print("\n Wartosci po zastosowaniu algorytmu to: \n Sx = ", np.matrix.round(Sx2g), "\n a =", np.around(a2,5))
polyn_write(m1,np.around(a1,5),m1)
polyn_write(m2,np.around(a2,5),m2)
print("\n Liczba punktów pomiarowych to:", np.size(X))
print("\n WIelkosci statystystyczne wyglądają następująco: \n RMSD1=", RMSD1, "\n RMSD2=", RMSD2,"\n E1=", E1, "\n E2=", E2)


plt.figure()
plt.plot(X, Y, color='#DB7093',linewidth=2, alpha=0.5, label='Dane pomiarowe')     
plt.plot(X,P1, color='blue', linewidth=2, alpha=0.5, label='Wielomian dopasowany m=5')
plt.plot(X,P2, color='red', linewidth=2, alpha=0.5, label='Wielomian dopasowany m=10')
plt.scatter(0, P1[0], marker='+', color='#8a134e',alpha=1.0, label='Punkt r1=0') 
plt.scatter(0, P2[0], marker='*', color='black',alpha=1.0, label='Punkt r2=0') 
#plt.ylim(17,10)
plt.xlabel('Odległosć [cm]/1e10')
plt.ylabel('Cisnienie [bar]/1e11')  
plt.legend(loc='upper right')
plt.title("Rozkład cisnienia Słońca w zależnosci od promienia liczone Gaussem")
plt.savefig("Rozklad_gauss.png")  
plt.show()


'''
Lepsze jest dopasowanie dla m=10, co wynika z wartosci RMSD oraz E
'''



