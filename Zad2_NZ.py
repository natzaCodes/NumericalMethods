"""
@author: Natalia Zalewska
Zadanie Zaliczeniowe 2 - wyliczanie wartosci zadanego wielomianu z użyciem schematu Hornera
"""


import numpy as np
import matplotlib.pyplot as plt


'''podanie stopnia wielomianu i współczynników'''

print("Podaj stopień wielomianu: ")
N = int(input())
while N<5:
    print("Stopień ma być większy lub równy 5: ")
    N = int(input())

n=0
atab = []

while n<=N:
    print("Podaj", n, "współczynnik dla wielomianu", N, "stopnia:")
    a = float(input())
    atab.append(a)
    n+=1

# przyda sie do wykresu
napis = "Wspolczynniki to: \n"
n=0
while n<=N:
    napis = napis + str(atab[n]) + ', '
    n+=1
    
    
    
    
''' Wypisanie wzoru'''
print("\n Twój wielomian to:")
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
 
    
 
''' Deklaracja dziedziny i obliczenia funkcji w dziedzinie'''

print("\n \n Podaj dziedzinę dla wykresu: [xmin; xmax]. Podaj xmin:")
xmin = np.float32(input())
print("\n Podaj xmax:")
xmax = np.float32(input())   
kroki = int((xmax - xmin)*100)

X=np.linspace(xmin, xmax, kroki)
Ptab = []

dane = open('dane.txt', "w")
dane.write("x\tP(x)\n")



''' + Zapis do pliku + wykres, schemat Hornera'''
for x in X:
    n=N 
    P=0
    while n>=0:
        P = P*x + atab[n]
        n-=1
    
    Ptab.append(P)
    dane.write(str(x))
    dane.write("\t")
    dane.write(str(P))
    dane.write("\n")

dane.close()
  
plt.scatter(X, Ptab, s=5)
plt.title("Wielomian stopnia {0}. \n {1} \n".format(N, napis),fontsize=10, fontweight=0) 
plt.xlabel('x')
plt.ylabel('Wartosc')  
plt.savefig("Wykres_NZ.png",bbox_inches='tight', pad_inches=0) 




''' Dzielenie przez dwumian'''

print("Podaj x0.")
x0 = np.float32(input())
n=N 
P=0
btab = np.empty([N+1])


while n>=0:
    if n==N:
        b=atab[n]
    else:
        b= atab[n] + btab[n+1]*x0
        
    P = P*x0 + atab[n]
    btab[n] = b
    n-=1    

print("\n Wynik dzielenia to: \n")
if x0>0:
    print("P(x) = (x -", x0, ') (', end="")
elif x0 == 0:
    print("P(x) = (x) (", end="")
else :
    print("P(x) = (x +", -x0, ') (', end="")

n=1
while n<=N:
      
    if btab[n]<0:
        print(" ", end="")
    elif btab[n] == 0 or n==1:
        print(end="") 
    else:
        print(' + ', end="")
        
        
    if btab[n] != 0 and btab[n] !=1 and btab[n] !=(-1) and n>=2:
        print(btab[n], end="")
    elif btab[n] !=0 and n==1:
        print(btab[n], end="")
    elif btab[n] ==1 and n>=2:
        print(end="")
    elif btab[n] == (-1) and n>=2:
        print("-", end="")
   
    if n == 1 or n==-1 or btab[n] == 0:
        print(end="")
    elif n == 2:
        print('x', end="")
    else:
        print('x^', end="") 
        print(n-1, end="")
        
    n+=1
        
print(")", end="")   
if btab[0]<0 or btab[0] == 0:
    print(" ", end="")
else:
    print(' + ', end="")   
    
print(btab[0])
    

''' Miłej zabawy :) '''
    

















