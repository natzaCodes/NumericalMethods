"""
Created on Fri Apr 22 14:22:13 2022
@author: Natalia
Zadanie zaliczeniowe 6
"""
import numpy as np
import sys


def Gauss_Jordan(macierz_glowna, data,czy_liczyc_macierz_odwrtona):
    
    matrix=np.copy(data)  
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
    
def Dolittle(data):
    N = len(data)
    L = np.identity(N)
    u = np.zeros((N,N))
    
    for j in range(N):
        for k in range(j,N):
            u[j,k]=data[j,k]
            for l in range(j):
                u[j,k]-=L[j,l]*u[l,k]
        for k in range(j+1,N):
            L[k,j]=data[k,j]
            for l in range(j):
                L[k,j]-=L[k,l]*u[l,j]
            if u[j,j]==0:
                sys.exit()
            L[k,j]/=u[j,j]
            
    return L,u
   
def LU(macierz_glowna,data,czy_liczyc_wyznacznik):
    L,U =Dolittle(data)
    b=data[:,-1]
    y=np.zeros(len(b))
    N = len(y)
    x=np.zeros(N)
    
    for j in range(N):
        if U[j,j]==0:
            sys.exit()
        y[j] = b[j]
        for k in range(j):
            y[j]-=L[j,k]*y[k]
    for j in range(N-1,-1,-1):   
           x[j] = y[j]
           for k in range(j+1,N):
               x[j]-=U[j,k]*x[k]

           x[j]/=U[j,j]
           
    if czy_liczyc_wyznacznik==True:
        wyznacznik = U.diagonal().prod() #*L.diagonal().prod()           
        return x, np.around(wyznacznik)
    
    return x
            
#data = np.array([[1.,2.,3.,4.],[4.,1.,2.,3.],[3.,1.,2.,4.]])

'''Main code'''

#dane z pliku
with open('uklad.txt') as plik:
    tab = [list(map(float, wiersz.split(' '))) for wiersz in plik]
uklad = np.array(tab)
macierz_glowna=uklad[:,0:len(uklad[1,:])-1] 


# wywołanie GJ 
odpowiedz_mo = input("Czy ma liczyć macierz odwrotną?")
if odpowiedz_mo == 'Tak' or odpowiedz_mo=='tak':
    GJ = Gauss_Jordan(macierz_glowna,uklad,czy_liczyc_macierz_odwrtona=True)
    print("\n \n Macierz główna: \n", np.around(GJ[0]),"\n \n x=",np.around(GJ[1]),"\n \n A^(-1): \n",np.around(GJ[2]))
else:    
    GJ = Gauss_Jordan(macierz_glowna,uklad,czy_liczyc_macierz_odwrtona=False)
    print("\n \n Macierz główna: \n", np.around(GJ[0]),"\n \n x=",np.around(GJ[1]))


# wywołanie LU
odpowiedz_w = input("Czy ma liczyć wyznacznik?")
if odpowiedz_w == 'Tak' or odpowiedz_w=='tak':
    LU = LU(macierz_glowna,uklad,czy_liczyc_wyznacznik=True)
    print("\n \n x=", np.around(LU[0]),"\n \n det A=",np.around(LU[1]))
else:    
    LU = LU(macierz_glowna,uklad,czy_liczyc_wyznacznik=True)
    print("\n \n x=", np.around(LU[0]))

      

       
            
            
            
            
            
            
            
            
            