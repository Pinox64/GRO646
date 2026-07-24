import Variables as var
import numpy as np
import math
import Ts_finder
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp


##############
# Paramètres #
##############
cas = 3
if cas == 1:
    temp_choice = 0 #50K
    wind = 0

elif cas == 2:
    temp_choice = 0 #50K
    wind = var.vent

elif cas == 3:
    temp_choice = 4 #250K
    wind = 0

elif cas == 4:
    temp_choice = 4 #250K
    wind = var.vent 


vec_actif = True
switch_actif = True

T_inf = var.T[temp_choice]

deltaT = abs(T_inf - var.TbPerm)

Lc = var.D

#first guess
Ts_first = var.TbPerm

r_out = var.D/2
r_in = var.Db/2

def Q (Ts, var):

    tf = (Ts + T_inf)/2

    Beta = 1/tf

    Ra_D = (var.g * Beta * (Ts - T_inf) * (Lc**3) * var.Pr[temp_choice]) / (var.nu[temp_choice]**2)

    if cas == 1 or cas == 3:
        Nu_sphere_Nat = 2 + (0.589*(Ra_D**(1/4)))/(1 + (0.469/var.Pr[temp_choice])**(9/16))**(4/9)
        h = (Nu_sphere_Nat * var.k[temp_choice])/Lc

    elif cas == 2 or cas == 4:

        Re=wind*var.D/var.nu[temp_choice]
        u_inf=var.rho[temp_choice]*wind*var.D/Re

        Close_index = np.argmin((np.abs(x - Ts) for x in var.T))
        Re_s=wind*var.D/var.nu[Close_index]
        u_s= var.rho[Close_index]*wind*var.D/Re_s

        Nuforce=2+(0.4*np.sqrt(Re)+0.06*(Re**(2/3)))*(var.Pr[temp_choice]**0.4)*((u_inf/u_s)**0.25)
        h=Nuforce*var.k[temp_choice]/var.D
    
    #####################################
    # Transferts de chaleurs mis en jeu #
    #####################################

    Qgaz=0.98
    
    # Puissance absorbée du rayonnement solar(W)
    Qsun = var.sigma * var.eps_s * (var.A_s / 2) * (var.r_s/(var.r_max if temp_choice == 0 else var.r_min))**2 * var.T_sun**4

    # Puissance rayonnée par le corps du robot (W)
    Qrad = var.eps_s * var.sigma * (var.A_s-var.A_vec) * (Ts**4 - T_inf**4)

    # Puissance rayonnée par le vec (w)
    Qrad_vec = (var.eps_vec if vec_actif else var.eps_s) * var.sigma * (var.A_vec) * (Ts**4 - T_inf**4)

    Qconv = h*var.A_s * (Ts-T_inf)
   
    Qconds = (4 * math.pi * var.k_in * r_in * r_out / (r_out - r_in)) * (var.TbPerm - Ts)

    # pas utilisé
    Qcondg=var.k_in*var.A_cont*(var.TbPerm -T_inf)/var.x_sw
    # pas utilisé
    Qcondsw= (var.k_sw*var.A_sw*(var.TbPerm-T_inf)/var.x_sw) if switch_actif else 0

    return Qsun, Qrad, Qrad_vec, Qconv, Qconds, Qcondg, Qcondsw, Qgaz
    

def Q_t(Ts, var):

    (
        Qsun,
        Qrad,
        Qrad_vec,
        Qconv,
        Qconds,
        Qcondg,
        Qcondsw,
        Qgaz
    ) = Q(Ts, var)

    Q_t = Qsun - Qrad - Qrad_vec - Qconv + Qconds
    return Q_t

Ts = Ts_finder.newton_raphson(Ts_first, lambda Ts: Q_t(Ts, var))

(
    Qsun,
    Qrad,
    Qrad_vec,
    Qconv,
    Qconds,
    Qcondg,
    Qcondsw,
    Qgaz
) = Q(Ts, var)

print(f"Ts: {Ts}")
print(f"Qsun: {Qsun}")
print(f"Qrad: {Qrad}")
print(f"Qrad_vec: {Qrad_vec}")
print(f"Qconv: {Qconv}")
print(f"Qconds: {Qconds}")
print(f"Qcondg: {Qcondg}")
print(f"Qcondsw: {Qcondsw}")
print(f"Qgaz: {Qgaz}")
print(f"Qt: {Q_t(Ts, var)}")

#Ein
print(f"Ein: {Qsun + Qconds}")

#Eout = Q_loss
print(f"Eout: {Qrad + Qconv + Qrad_vec}")
