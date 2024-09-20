"""
@author: Natalia
Date: 3.05.2022
Least Square Method
"""

import numpy as np
import math
import matplotlib.pyplot as plt

def dopasowanie(m_w, log_P, sigmaw):
# szukam mwdop = a log10P + b
    S=0
    Sx=0
    Sy=0
    Sxx=0
    Syy=0
    Sxy=0
   
    for i in range(len(m_w)):               
        S += 1/sigmaw**2 
        Sx += log_P[i]/sigmaw**2
        Sy += m_w[i]/sigmaw**2
        Sxx += log_P[i]**2/sigmaw**2
        Syy += m_w[i]**2/sigmaw**2
        Sxy += log_P[i]*m_w[i]/sigmaw**2
        
    delta = S*Sxx - (Sx)**2
    r = (S*Sxy-Sx*Sy)/np.sqrt((S*Sxx-Sx**2)*(S*Syy-Sy**2))
    
    alpha = (S*Sxy-Sx*Sy)/delta
    beta = (Sxx*Sy-Sx*Sxy)/delta
    
    sigmaa = np.sqrt(S/delta)
    sigmab = np.sqrt(Sxx/delta)
    
    #wartosci funkcji dopasowanej
    f=0
    f_tab = []
    Chi2=0
    RMSD=0
    for i in range(len(m_w)): 
        f =  alpha*log_P[i] + beta
        Chi2 += (m_w[i]-f)**2/sigmaw**2
        RMSD += (m_w[i]-f)**2
        f_tab.append(f)
    RMSD = np.sqrt(RMSD/len(m_w))
    Chi2=Chi2/len(m_w)
    
    return (f_tab, alpha, beta, sigmaa, sigmab, Chi2, RMSD, r)  

'''dane z pliku'''

cepF = open("cep.txt", "r")
linescepF = cepF.readlines()
cepF.close()

m_I = []
m_V = []
Period = []

#podział na listy
for i in linescepF:   
    temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11 = i.strip().split(" ")
   
    m_I.append(str(temp2))
    m_V.append(str(temp3))
    Period.append(str(temp4))

#wywalanie "-" i zmiana str na float
i=0
while i < len(m_I):
    if m_I[i] == "-" or m_V[i] == "-" or Period[i] == "-":
        del m_I[i]
        del m_V[i]
        del Period[i]
    else:
        m_I[i] = float(m_I[i])
        m_V[i] = float(m_V[i])
        Period[i] = float(Period[i])
        i+=1 


'''indeks Wesenheit i log(P)'''

sigmaI = 0.02 #mag
sigmaV = 0.02 #mag
sigmaw = np.sqrt((2.55**2*sigmaI**2) + (1.55**2*sigmaV**2)) #z reguły propagacji bledów
m_w = []
log_P = []

for i in range(len(m_I)):
       mw = m_I[i] - 1.55*(m_V[i]-m_I[i])
       logP = math.log10(Period[i])
       m_w.append(mw)
       log_P.append(logP)


''' Oddalone punkty od dopasowania '''

f_tab, alpha, beta, sigmaa, sigmab, Chi2, rmsd, r = dopasowanie(m_w, log_P, sigmaw)
print("Liczba punktów uwzględnianych:", len(m_w))
print("\n współczynnik alpha:", round(alpha, 4),"+-", round(sigmaa, 4))
print("\n współczynnik beta:", round(beta, 4),"+-", round(sigmab,4))
print("\n Wartosc Chi^2/dof:", round(Chi2, 4))
print("\n Wartosc RMSD:", round(rmsd,4))
print("\n Współczynnik korelacji liniowej Persona:", round(r, 4))


N_out = 1
i=0
m_discards =[]
P_discards = []

while N_out>0:
    N_out =0
    n=0
    i=0
    while i < len(m_w):
        roznica = np.abs(m_w[i]-f_tab[i])
        if roznica > 3*rmsd:
            m_discards.append(m_w[i])
            P_discards.append(log_P[i])
            del m_w[i]
            del log_P[i]
            del f_tab[i]
            n+=1
            N_out+=1
        else:
            i+=1
            
    f_tab, alpha, beta, sigmaa, sigmab, Chi2, rmsd, r = dopasowanie(m_w, log_P, sigmaw)
    print("\n \nLiczba punktów uwzględnianych:", i)
    print("\n współczynnik alpha:", round(alpha, 4),"+-", round(sigmaa, 4))
    print("\n współczynnik beta:", round(beta, 4),"+-", round(sigmab,4))
    print("\n Wartosc Chi^2/dof:", round(Chi2, 4))
    print("\n Wartosc RMSD:", round(rmsd,4))
    print("\n Współczynnik korelacji liniowej Persona:", round(r, 4))
            
plt.scatter(log_P, m_w, marker='+', color='#8a134e',alpha=1.0, label='jasnosc Wesenheit #choices') 
plt.scatter(P_discards, m_discards, marker='*', color='#DB7093',alpha=0.5, label='jasnosc Wesenheit #discards')     
plt.plot(log_P, f_tab, color='blue', linewidth=1, alpha=0.5, label='prosta dopasowana')
plt.ylim(17,10)
plt.xlabel('Period [day] in log10 scale')
plt.ylabel('Jasnosc Wesenheit [mag]')  
plt.legend(loc='lower right')
plt.title("Jasnosć Wesenheit dla cefeid")
plt.savefig("Wesenheit.png")  


'''
Wraz z iteracjami sigma alfpha i sigma beta nieznacznie się zwiększają, Zredukowane Chi2 dąży do 1 i przy ostatniej iteracji ~1, 
co oznacza, że dopasowanie jest coraz dokładniejsze a na koniec jest prawie idealne.
błąd sredniokwadratowy maleje co jest dobrym sygnałem, gdyż rozrzut danych wokół dopasowania maleje,
a współczynnik korelacji liniowej dąży do -1, z czego można wywnioskować, że zachodzi korelacja negatywna 
'''


