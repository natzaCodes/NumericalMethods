"""
Created on Fri May 20 15:29:28 2022
@author: Natalia
Zadanie zaliczeniowe 9
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

def trapezy(f,a,b,eps):
    h0=b-a
    I0=h0*(f(a)+f(b))/2
    h1=h0/2
    I1=I0/2 + h1*f((a+b)/2)
    n=4
    while(abs(I1-I0)/abs((I0)+1e-15)>eps):
        I0=I1
        h1=h1/2
        I1=I1/2
        for i in range (1,n,2):
            I1=I1+h1*(f(a+i*h1))
        n=n*2
    return I1

def Simpson(f,a,b,eps):
    
    h0=(b-a)/2
    I0=h0*(f(a)+4*f(a+h0)+f(b))/3
    h1=h0/2
    I1=I0/2+2*h1*(2*(f(a+h1)+f(a+3*h1))-f(a+2*h1))/3
    
    n=8
    
    while(abs(I0-I1)/(abs(I0)+1e-15)>eps):
        I0=I1
        h1/=2
        I1/=2
        np, p= 0, 0 #nieparzyste, parzyste
        
        for i in range(1, n, 2):
            np+=f(a+i*h1)
        for i in range(2, n-1, 4):
            p+=f(a+i*h1)
        I1+=h1*(4*np-2*p)/3
        n*=2
        
    return I1
    
    n=8
    #dokończy pisac
    
    return I1

def phi(M):
    #calkowita liczba gwiazd w probce
    A1=0.256274
    A2=0.137334
    A3=0.137334

    if M>=0.08 and M<0.5:
        return A1*M**(-1.3)
    elif M>=0.5 and M<1.0:
        return A2*M**(-2.2)
    elif M>=1.0 and M<=150:
        return A3*M**(-2.7)

def IMF(M):
    return M*phi(M)

def lumin(M):   
    if M>=0.08 and M<1.0:
        return phi(M)*M**(5)
    if M>=1.0 and M<=150.0:
        return phi(M)*M**(3)

def rad(M):
    if M>=0.08 and M<1.5:
        return phi(M)*M**(0.4)
    if M>=1.5 and M<=150.0:
        return phi(M)*M**(0.8)
    

''' Main code'''
N=1e11 #calkowita liczba gwiazd w Drodze Mlecznej
eps=1e-5


''' 1 przedział'''
Mpa=0.08
Mka=150.0
Msimpa=Simpson(IMF,Mpa,Mka,eps)/Simpson(phi,Mpa,Mka,eps)
Mgla=GL(IMF,Mpa,Mka)/GL(phi,Mpa,Mka)
Mglra=GL_rec(IMF,Mpa,Mka,eps,GL(IMF,Mpa,Mka))/GL_rec(phi,Mpa,Mka,eps,GL(phi,Mpa,Mka))
Mtrapa= trapezy(IMF,Mpa,Mka,eps)/trapezy(phi,Mpa,Mka,eps)

Nsimpa=Simpson(phi,Mpa,Mka,eps)/Simpson(phi,Mpa,Mka,eps)*N
Ngla=GL(phi,Mpa,Mka)/GL(phi,Mpa,Mka)*N
Nglra=GL_rec(phi,Mpa,Mka,eps,GL(phi,Mpa,Mka))/GL_rec(phi,Mpa,Mka,eps,GL(phi,Mpa,Mka))*N
Ntrapa= trapezy(phi,Mpa,Mka,eps)/trapezy(phi,Mpa,Mka,eps)*N

Lsimpa=Simpson(lumin,Mpa,Mka,eps)/Simpson(phi,Mpa,Mka,eps)
Lgla=GL(lumin,Mpa,Mka)/GL(phi,Mpa,Mka)
Lglra=GL_rec(lumin,Mpa,Mka,eps,GL(lumin,Mpa,Mka))/GL_rec(phi,Mpa,Mka,eps,GL(phi,Mpa,Mka))
Ltrapa= trapezy(lumin,Mpa,Mka,eps)/trapezy(phi,Mpa,Mka,eps)

Rsimpa=Simpson(rad,Mpa,Mka,eps)/Simpson(phi,Mpa,Mka,eps)
Rgla=GL(rad,Mpa,Mka)/GL(phi,Mpa,Mka)
Rglra=GL_rec(rad,Mpa,Mka,eps,GL(rad,Mpa,Mka))/GL_rec(phi,Mpa,Mka,eps,GL(phi,Mpa,Mka))
Rtrapa= trapezy(rad,Mpa,Mka,eps)/trapezy(phi,Mpa,Mka,eps)


''' 2 przedział '''
Mpb=0.08
Mkb=1.0
Msimpb=Simpson(IMF,Mpb,Mkb,eps)/Simpson(phi,Mpb,Mkb,eps)
Mglb=GL(IMF,Mpb,Mkb)/GL(phi,Mpb,Mkb)
Mglrb=GL_rec(IMF,Mpb,Mkb,eps,GL(IMF,Mpb,Mkb))/GL_rec(phi,Mpb,Mkb,eps,GL(phi,Mpb,Mkb))
Mtrapb= trapezy(IMF,Mpb,Mkb,eps)/trapezy(phi,Mpb,Mkb,eps)

Nsimpb=Simpson(phi,Mpb,Mkb,eps)/Simpson(phi,Mpa,Mka,eps)*N
Nglb=GL(phi,Mpb,Mkb)/GL(phi,Mpa,Mka)*N
Nglrb=GL_rec(phi,Mpb,Mkb,eps,GL(phi,Mpb,Mkb))/GL_rec(phi,Mpa,Mka,eps,GL(phi,Mpa,Mka))*N
Ntrapb= trapezy(phi,Mpb,Mkb,eps)/trapezy(phi,Mpa,Mka,eps)*N

Lsimpb=Simpson(lumin,Mpb,Mkb,eps)/Simpson(phi,Mpb,Mkb,eps)
Lglb=GL(lumin,Mpb,Mkb)/GL(phi,Mpb,Mkb)
Lglrb=GL_rec(lumin,Mpb,Mkb,eps,GL(lumin,Mpb,Mkb))/GL_rec(phi,Mpb,Mkb,eps,GL(phi,Mpb,Mkb))
Ltrapb= trapezy(lumin,Mpb,Mkb,eps)/trapezy(phi,Mpb,Mkb,eps)

Rsimpb=Simpson(rad,Mpb,Mkb,eps)/Simpson(phi,Mpb,Mkb,eps)
Rglb=GL(rad,Mpb,Mkb)/GL(phi,Mpb,Mkb)
Rglrb=GL_rec(rad,Mpb,Mkb,eps,GL(rad,Mpb,Mkb))/GL_rec(phi,Mpb,Mkb,eps,GL(phi,Mpb,Mkb))
Rtrapb= trapezy(rad,Mpb,Mkb,eps)/trapezy(phi,Mpb,Mkb,eps)


''' 3 przedział '''
Mpc=1.0
Mkc=8.0
Msimpc=Simpson(IMF,Mpc,Mkc,eps)/Simpson(phi,Mpc,Mkc,eps)
Mglc=GL(IMF,Mpc,Mkc)/GL(phi,Mpc,Mkc)
Mglrc=GL_rec(IMF,Mpc,Mkc,eps,GL(IMF,Mpc,Mkc))/GL_rec(phi,Mpc,Mkc,eps,GL(phi,Mpc,Mkc))
Mtrapc= trapezy(IMF,Mpc,Mkc,eps)/trapezy(phi,Mpc,Mkc,eps)

Nsimpc=Simpson(phi,Mpc,Mkc,eps)/Simpson(phi,Mpa,Mka,eps)*N
Nglc=GL(phi,Mpc,Mkc)/GL(phi,Mpa,Mka)*N
Nglrc=GL_rec(phi,Mpc,Mkc,eps,GL(phi,Mpb,Mkb))/GL_rec(phi,Mpa,Mka,eps,GL(phi,Mpa,Mka))*N
Ntrapc= trapezy(phi,Mpc,Mkc,eps)/trapezy(phi,Mpa,Mka,eps)*N

Lsimpc=Simpson(lumin,Mpc,Mkc,eps)/Simpson(phi,Mpc,Mkc,eps)
Lglc=GL(lumin,Mpc,Mkc)/GL(phi,Mpc,Mkc)
Lglrc=GL_rec(lumin,Mpc,Mkc,eps,GL(lumin,Mpc,Mkc))/GL_rec(phi,Mpc,Mkc,eps,GL(phi,Mpc,Mkc))
Ltrapc= trapezy(lumin,Mpc,Mkc,eps)/trapezy(phi,Mpc,Mkc,eps)

Rsimpc=Simpson(rad,Mpc,Mkc,eps)/Simpson(phi,Mpc,Mkc,eps)
Rglc=GL(rad,Mpc,Mkc)/GL(phi,Mpc,Mkc)
Rglrc=GL_rec(rad,Mpc,Mkc,eps,GL(rad,Mpc,Mkc))/GL_rec(phi,Mpc,Mkc,eps,GL(phi,Mpc,Mkc))
Rtrapc= trapezy(rad,Mpc,Mkc,eps)/trapezy(phi,Mpc,Mkc,eps)


''' 4 przedział '''
Mpd=8.0
Mkd=21.0
Msimpd=Simpson(IMF,Mpd,Mkd,eps)/Simpson(phi,Mpd,Mkd,eps)
Mgld=GL(IMF,Mpd,Mkd)/GL(phi,Mpd,Mkd)
Mglrd=GL_rec(IMF,Mpd,Mkd,eps,GL(IMF,Mpd,Mkd))/GL_rec(phi,Mpd,Mkd,eps,GL(phi,Mpd,Mkd))
Mtrapd= trapezy(IMF,Mpd,Mkd,eps)/trapezy(phi,Mpd,Mkd,eps)

Nsimpd=Simpson(phi,Mpd,Mkd,eps)/Simpson(phi,Mpa,Mka,eps)*N
Ngld=GL(phi,Mpd,Mkd)/GL(phi,Mpa,Mka)*N
Nglrd=GL_rec(phi,Mpd,Mkd,eps,GL(phi,Mpd,Mkd))/GL_rec(phi,Mpa,Mka,eps,GL(phi,Mpa,Mka))*N
Ntrapd= trapezy(phi,Mpd,Mkd,eps)/trapezy(phi,Mpa,Mka,eps)*N

Lsimpd=Simpson(lumin,Mpd,Mkd,eps)/Simpson(phi,Mpd,Mkd,eps)
Lgld=GL(lumin,Mpd,Mkd)/GL(phi,Mpd,Mkd)
Lglrd=GL_rec(lumin,Mpd,Mkd,eps,GL(lumin,Mpd,Mkd))/GL_rec(phi,Mpd,Mkd,eps,GL(phi,Mpd,Mkd))
Ltrapd= trapezy(lumin,Mpd,Mkd,eps)/trapezy(phi,Mpd,Mkd,eps)

Rsimpd=Simpson(rad,Mpd,Mkd,eps)/Simpson(phi,Mpd,Mkd,eps)
Rgld=GL(rad,Mpd,Mkd)/GL(phi,Mpd,Mkd)
Rglrd=GL_rec(rad,Mpd,Mkd,eps,GL(rad,Mpd,Mkd))/GL_rec(phi,Mpd,Mkd,eps,GL(phi,Mpd,Mkd))
Rtrapd= trapezy(rad,Mpd,Mkd,eps)/trapezy(phi,Mpd,Mkd,eps)


''' 5 przedział '''
Mpe=21.0
Mke=150.0
Msimpe=Simpson(IMF,Mpe,Mke,eps)/Simpson(phi,Mpe,Mke,eps)
Mgle=GL(IMF,Mpe,Mke)/GL(phi,Mpe,Mke)
Mglre=GL_rec(IMF,Mpe,Mke,eps,GL(IMF,Mpe,Mke))/GL_rec(phi,Mpe,Mke,eps,GL(phi,Mpe,Mke))
Mtrape= trapezy(IMF,Mpe,Mke,eps)/trapezy(phi,Mpe,Mke,eps)

Nsimpe=Simpson(phi,Mpe,Mke,eps)/Simpson(phi,Mpa,Mka,eps)*N
Ngle=GL(phi,Mpe,Mke)/GL(phi,Mpa,Mka)*N
Nglre=GL_rec(phi,Mpe,Mke,eps,GL(phi,Mpe,Mke))/GL_rec(phi,Mpa,Mka,eps,GL(phi,Mpa,Mka))*N
Ntrape= trapezy(phi,Mpe,Mke,eps)/trapezy(phi,Mpa,Mka,eps)*N

Lsimpe=Simpson(lumin,Mpe,Mke,eps)/Simpson(phi,Mpe,Mke,eps)
Lgle=GL(lumin,Mpe,Mke)/GL(phi,Mpe,Mke)
Lglre=GL_rec(lumin,Mpe,Mke,eps,GL(lumin,Mpe,Mke))/GL_rec(phi,Mpe,Mke,eps,GL(phi,Mpe,Mke))
Ltrape= trapezy(lumin,Mpe,Mke,eps)/trapezy(phi,Mpe,Mke,eps)

Rsimpe=Simpson(rad,Mpe,Mke,eps)/Simpson(phi,Mpe,Mke,eps)
Rgle=GL(rad,Mpe,Mke)/GL(phi,Mpe,Mke)
Rglre=GL_rec(rad,Mpe,Mke,eps,GL(rad,Mpe,Mke))/GL_rec(phi,Mpe,Mke,eps,GL(phi,Mpe,Mke))
Rtrape= trapezy(rad,Mpe,Mke,eps)/trapezy(phi,Mpe,Mke,eps)



zaokrl=5
print(" \n \n Srednia masa gwiazd:")
print("przedziały: \t [0.08;150]Ms \t [0.08;1.0]Ms \t [1;8]Ms \t [8;21]Ms \t [21;150]Ms")
print("met. G-L podstawowa: \t",np.round(Mgla,zaokrl),"\t",np.round(Mglb,zaokrl),"\t ",np.round(Mglc,zaokrl),"\t",np.round(Mgld,zaokrl),"\t",np.round(Mgle,zaokrl))
print("met. G-L iteracyjna: \t",np.round(Mglra,zaokrl),"\t",np.round(Mglrb,zaokrl),"\t",np.round(Mglrc,zaokrl),"\t", np.round(Mglrd,zaokrl),"\t",np.round(Mglre,zaokrl))
print("met. trapezów: \t",np.round(Mtrapa,zaokrl),"\t",np.round(Mtrapb,zaokrl),"\t",np.round(Mtrapc,zaokrl),"\t",np.round(Mtrapd,zaokrl),"\t",np.round(Mtrape,zaokrl))
print("met. Simpsona:\t",np.round(Msimpa,zaokrl),"\t",np.round(Msimpb,zaokrl),"\t",np.round(Msimpc,zaokrl),"\t", np.round(Msimpd,zaokrl)," \t",np.round(Msimpe,zaokrl))

zaokrl=1
print("\n \n Liczba gwiazd:")
print("przedziały: \t [0.08;150]Ms \t [0.08;1.0]Ms \t [1;8]Ms \t [8;21]Ms \t [21;150]Ms")
print("met. G-L podstawowa: \t",np.round(Ngla,zaokrl),"\t",np.round(Nglb,zaokrl),"\t ",np.round(Nglc,zaokrl),"\t",np.round(Ngld,zaokrl),"\t",np.round(Ngle,zaokrl))
print("met. G-L iteracyjna: \t",np.round(Nglra,zaokrl),"\t",np.round(Nglrb,zaokrl),"\t",np.round(Nglrc,zaokrl),"\t", np.round(Nglrd,zaokrl),"\t",np.round(Nglre,zaokrl))
print("met. trapezów: \t",np.round(Ntrapa,zaokrl),"\t",np.round(Ntrapb,zaokrl),"\t",np.round(Ntrapc,zaokrl),"\t",np.round(Ntrapd,zaokrl),"\t",np.round(Ntrape,zaokrl))
print("met. Simpsona:\t",np.round(Nsimpa,zaokrl),"\t",np.round(Nsimpb,zaokrl),"\t",np.round(Nsimpc,zaokrl),"\t", np.round(Nsimpd,zaokrl)," \t",np.round(Nsimpe,zaokrl))

zaokrl=5
print("\n \n Średnia jasnosć gwiazdy:")
print("przedziały: \t [0.08;150]Ms \t [0.08;1.0]Ms \t [1;8]Ms \t [8;21]Ms \t [21;150]Ms")
print("met. G-L podstawowa: \t",np.round(Lgla,zaokrl),"\t",np.round(Lglb,zaokrl),"\t ",np.round(Lglc,zaokrl),"\t",np.round(Lgld,zaokrl),"\t",np.round(Lgle,zaokrl))
print("met. G-L iteracyjna: \t",np.round(Lglra,zaokrl),"\t",np.round(Lglrb,zaokrl),"\t",np.round(Lglrc,zaokrl),"\t", np.round(Lglrd,zaokrl),"\t",np.round(Lglre,zaokrl))
print("met. trapezów: \t",np.round(Ltrapa,zaokrl),"\t",np.round(Ltrapb,zaokrl),"\t",np.round(Ltrapc,zaokrl),"\t",np.round(Ltrapd,zaokrl),"\t",np.round(Ltrape,zaokrl))
print("met. Simpsona:\t",np.round(Lsimpa,zaokrl),"\t",np.round(Lsimpb,zaokrl),"\t",np.round(Lsimpc,zaokrl),"\t", np.round(Lsimpd,zaokrl)," \t",np.round(Lsimpe,zaokrl))

zaokrl=5
print("\n \n Średni promień gwiazdy:")
print("przedziały: \t [0.08;150]Ms \t [0.08;1.0]Ms \t [1;8]Ms \t [8;21]Ms \t [21;150]Ms")
print("met. G-L podstawowa: \t",np.round(Rgla,zaokrl),"\t",np.round(Rglb,zaokrl),"\t ",np.round(Rglc,zaokrl),"\t",np.round(Rgld,zaokrl),"\t",np.round(Rgle,zaokrl))
print("met. G-L iteracyjna: \t",np.round(Rglra,zaokrl),"\t",np.round(Rglrb,zaokrl),"\t",np.round(Rglrc,zaokrl),"\t", np.round(Rglrd,zaokrl),"\t",np.round(Rglre,zaokrl))
print("met. trapezów: \t",np.round(Rtrapa,zaokrl),"\t",np.round(Rtrapb,zaokrl),"\t",np.round(Rtrapc,zaokrl),"\t",np.round(Rtrapd,zaokrl),"\t",np.round(Rtrape,zaokrl))
print("met. Simpsona:\t",np.round(Rsimpa,zaokrl),"\t",np.round(Rsimpb,zaokrl),"\t",np.round(Rsimpc,zaokrl),"\t", np.round(Rsimpd,zaokrl)," \t",np.round(Rsimpe,zaokrl))


'''
#do testów
print(masa(1.0,8.0,eps,trapez))
print(masa(8.0,21.0,eps,trapez))
print(masa(21.0,150.0,eps,trapez))
print(masa(1.0,8.0,eps,trapez))
print(masa(8.0,21.0,eps,trapez))
print(masa(21.0,150.0,eps,trapez))

I0=GL(f,-1,1)
I=GL_rec(f,-1,1,1e-5,I0)
print(I0)
print(I)
'''





























