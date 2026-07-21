# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 12:05:58 2024

@author: lebm1401
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import style
style.use('bmh')

def Calcul_Tsp(Tsp,h):
    Tsc=85+273 #K
    Tinf=298.   #K
    Tsurr=298.  #K
    
    Rcontprimeprime=2.75e-4 #m^2.K/W
    Ac=2e-4  #m^2
    Rcont=Rcontprimeprime/Ac #K/W
    
    k=240  #W/m.K
    L=6e-3 #m
    W=20e-3 #m
    A=W*W #m^2
    Rcond=L/(k*A)  #K/W
    
    Rconv=1/(h*A)  #K/W
    
    eps=0.9
    sigma=5.67e-8 #W/m^2.K^4
    hr=eps*sigma*(Tsp+Tsurr)*(Tsp**2+Tsurr**2) #W/m^2.K
    Rray=1/(hr*A) #K/W
    
    Req=Rconv*Rray/(Rconv+Rray) #K/W
    
    Rtot=Rcont+Rcond+Req # K/W
    
    Q1=(Tsc-Tinf)/Rtot #W
    
    Q2=(Tsc-Tsp)/(Rcont+Rcond)  #W

       
    return (Q1-Q2)


def rBissection(func,a,b,precision,h):
    i=0
    fa=func(a,h)
    fb=func(b,h)
    while (b-a) > precision:
        i=i+1
        m=(a+b)/2
        fm=func(m,h)
        if (fm*fb) < 0:
            a=m
            fa=func(a,h)
        else:
            b=m
            fb=func(b,h)
    x=m
    return x


precision=0.0001
a=298. #K
b=85.+273 #K

h=4. #W/m^2.K


(Tsp)=rBissection(Calcul_Tsp,a,b,precision,h)
skel="\nPour h={} W/m^2.K, Tsp = {} K "
print(skel.format(h,Tsp))



vec_Tsp=np.linspace(a,b,1000)
vec_calcul_Tsp=Calcul_Tsp(vec_Tsp,h)

plt.figure()
plt.figure(figsize=(8,6))
plt.plot(vec_Tsp,vec_calcul_Tsp,'-o')
plt.xlabel('Tsp (K)')
plt.ylabel('[(Tsc-Tinf)/Rtot]-[(Tsc-Tsp)/(Rcont+Rcond)] (W)')
plt.grid(True)

