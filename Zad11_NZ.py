"""
Created on Thu Jun  9 23:27:51 2022
@author: Natalia
Zadanie domowe 11
"""
import numpy as np
import matplotlib.pyplot as plt


def f(x,y,z):
    return -(2/x)*z-y**3

def Rk4(f,h,t,x0,y0,z0):
    n = len(t)
    y = np.zeros(n)
    z_t = np.zeros(n)
    for i in range(n):
        k1s = f(x0, y0, z0)
        k2s = f((x0+h/2), (y0+z0*h/2), (z0+h/2*k1s))
        k3s = f((x0+h/2), (y0+(z0+h/2*k1s)*h/2), (z0+h/2*k2s))
        k4s = f((x0+h), (y0+h*(z0+h/2*k2s)), (z0+h*k3s))
    
        y[i] = y0 + h*z0+ h*h*(k1s+k2s+k3s)/6
        y0 = y[i]
        x0 += h
        z = z0 + h/6*(k1s+2*k2s+2*k3s+k4s)
        z_t[i]=z
        z0=z
    return y, z_t

def AB4(f,h,t,Rk4,tr,x0,y0,z0):
    n = len(t)
    y1,z1=Rk4(f,h,tr,x0,y0,z0)

    y = np.zeros(n)
    y[0]=y0
    y[1]=y1[0]
    y[2]=y1[1]
    y[3]=y1[2]

    z = np.zeros(n)
    z[0]=z0
    z[1]=z1[0]
    z[2]=z1[1]
    z[3]=z1[2]

    x = x0 + 4*h
    for i in range(4,n):
        
        #print(f(t[i-1],y[i-1],z[i-1]),f(t[i-2],y[i-2], z[i-2]), f(t[i-3],y[i-3],z[i-3]), f(t[i-4],y[i-4],z[i-4]),"\n")
        y[i] = y[i-1]+ h/24*( 55*z[i-1] -59*z[i-2]+ 37*z[i-3] -9*z[i-4] )
        z[i] = z[i-1]+ h/24*( 55*f((x-h),y[i-1],z[i-1]) -59*f((x-2*h),y[i-2],z[i-2])+ 37*f((x-3*h),y[i-3],z[i-3]) -9*f((x-4*h),y[i-4],z[i-4]) )
        y0 = y[i]
        x+=h
        
    return y

def PC(f,h,t,x0,y0,z0):
    
    n = len(t)
    y = np.zeros(n)
    
    for i in range(n):

        #przewidywania
        z = z0 + h*f(x0,y0,z0)
        y[i] = y0 + h*z0
        #poprawianie
        x = x0+h
        z= z0 + h/2*(f(x0,y0,z0)+ f(x,y[i],z))
        y[i] = y0+ h/2*(z0+z)
        
        y0 = y[i]
        x0 = x
        z0=z    
    
    return y

def PECE(f,h,t,Rk4,tRk4,x0,y0,z0):
    
    n = len(t)
    y1,z1=Rk4(f,h,tRk4,x0,y0,z0)

    y = np.zeros(n)
    y[0]=y0
    y[1]=y1[0]
    y[2]=y1[1]
    y[3]=y1[2]

    z = np.zeros(n)
    z[0]=z0
    z[1]=z1[0]
    z[2]=z1[1]
    z[3]=z1[2]

    x = x0 + 4*h

    for i in range(4,n):
        #przewidywania
        y[i] = y[i-1]+ h/24*( 55*z[i-1] -59*z[i-2]+ 37*z[i-3] -9*z[i-4] )
        z[i] = z[i-1]+ h/24*( 55*f((x-h),y[i-1],z[i-1]) -59*f((x-2*h),y[i-2],z[i-2])+ 37*f((x-3*h),y[i-3],z[i-3]) -9*f((x-4*h),y[i-4],z[i-4]) )
        y0 = y[i]
        
        #evaluate
        fz=f(x,y[i],z[i])
        
        #poprawa
        y[i] = y[i-1]+ h/24*( 9*fz +19*z[i-1]- 5*z[i-2] +z[i-3] )
        z[i] = z[i-1]+ h/24*( 9*fz +19*f((x-h),y[i-1],z[i-1])- 5*f((x-2*h),y[i-2],z[i-2]) +f((x-3*h),y[i-3],z[i-3]) )
        
        #evaluate
        fz=f(x,y[i],z[i])
        
        x+=h
        
    return y


h=1e-4
zakresRk4=np.arange(0.,4*h,h)

zakres=np.arange(0.,10.,h)
y00=np.zeros(len(zakres))
#y_poczatek, z_poczatek = Rk4(f,h,zakresRk4,1e-20,1.,0.)

yAB = AB4(f,h,zakres,Rk4,zakresRk4,1e-16,1.,0.)
yPC = PC(f,h,zakres,1e-16,1.,0.)
yPECE = PECE(f,h,zakres,Rk4,zakresRk4,1e-16,1.,0.)


dane = open('dane.txt', "w")
dane.write("Nr\tx\tyAB\tyPC\tyPECE\n")


''' + Zapis do pliku'''
zaokrl=5
for i in range(len(zakres-1)):

    dane.write(str(i))
    dane.write("\t")
    dane.write(str(np.round(zakres[i],zaokrl)))
    dane.write("\t")
    dane.write(str(np.round(yAB[i],zaokrl)))
    dane.write("\t")
    dane.write(str(np.round(yPC[i],zaokrl)))
    dane.write("\t")
    dane.write(str(np.round(yPECE[i],zaokrl)))
    dane.write("\n")

dane.close()

plt.figure()
plt.plot(zakres,yAB,'purple',linewidth=1, alpha=0.5, label='Adams-Bashforth IV rzędu')
plt.plot(zakres,yPC,'red',linewidth=1, alpha=0.5, label='PC')
plt.plot(zakres,yPECE,'blue',linewidth=1, alpha=0.5, label='PECE')
plt.plot(zakres,y00,'black',linewidth=1, alpha=1.0, label='$\Theta$=0')
plt.xlabel('$\Xi$')
plt.ylabel('$\Theta$')  
plt.legend(loc='upper right')
plt.title("Wykres rozwiazanego równania dla indeksu politropy 3")
plt.savefig("Full.png")  
plt.show()

zakresa=68958
zakresb=68978
plt.figure()
plt.plot(zakres[zakresa:zakresb],yAB[zakresa:zakresb],'purple',linewidth=1, alpha=0.5, label='Adams-Bashforth IV rzędu')
plt.plot(zakres[zakresa:zakresb],yPC[zakresa:zakresb],'red',linewidth=1, alpha=0.5, label='PC')
#plt.plot(zakres[zakresa:zakresb],yPECE[zakresa:zakresb],'blue',linewidth=1, alpha=0.5, label='PECE')
plt.plot(zakres[zakresa:zakresb],y00[zakresa:zakresb],'black',linewidth=2, alpha=0.5, label='$\Theta$=0')
plt.xlabel('$\Xi$')
plt.ylabel('$\Theta$')  
plt.legend(loc='upper right')
plt.title("Wykres rozwiazanego równania dla indeksu politropy 3, przecięcie z y=0 \n",fontsize=10, fontweight=0)
plt.savefig("Przyblizenie.png",bbox_inches='tight', pad_inches=0.1)  
plt.show()




















