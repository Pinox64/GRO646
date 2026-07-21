import Variables as var
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

temp_choice = 0
h_nat = 0

# Constant
Qsun = var.sigma * var.eps_s * (var.A_s / 2) * (var.r_s/(var.r_max if temp_choice == 0 else var.r_min))**2 * var.T_sun**4
Qgaz = 0.98

# Fonction de Ts
Qrad = lambda Ts: var.eps_s * var.sigma * (var.A_s-var.A_vec) * (Ts**4 - var.T[temp_choice]**4)
Qrad_vec = lambda Ts: var.eps_vec * var.sigma * (var.A_vec) * (Ts**4 - var.T[temp_choice]**4)
Qconv = lambda Ts, h_nat: h_nat*var.A_s * (Ts - var.T[temp_choice])
Qcondg = lambda Ts: var.k_in*var.A_cont*(var.T[temp_choice]-Ts)/var.x_ins

# Fonction de Tb
Qcondsw = lambda Tb: var.k_sw*var.A_sw*(Tb - var.T[temp_choice])/var.x_sw

# Fonction de Ts et Tb
Qconds = lambda Ts, Tb: (var.A_s*var.k_in/var.x_ins)*(Tb - Ts)

def bilan_surface(Ts_guess, Tb_actuelle):
    Ts = Ts_guess[0] # fsolve passe une liste
    
    chaleur_venant_du_corps = Qconds(Ts, Tb_actuelle)
    
    tf = (Ts + var.T[temp_choice])/2

    Beta = 1/tf

    Ra_D = (var.g * Beta * abs(Ts - var.T[temp_choice]) * (var.D**3) * var.Pr[temp_choice]) / var.nu[temp_choice]

    Nu_sphere_Nat = 2 + (0.589*(Ra_D**(1/6)))/(1 + (0.469/var.Pr[temp_choice])**(9/16))**(4/9)

    h_nat = (Nu_sphere_Nat * var.k[temp_choice])/var.D

    # Bilan = 0
    residu = chaleur_venant_du_corps + Qsun + Qgaz - Qrad(Ts) - Qconv(Ts, h_nat)
    return residu

def dTdt(t, Tb_actuelle):
    
    Ts_sol = fsolve(bilan_surface, x0=[Tb_actuelle], args=(Tb_actuelle,))
    Ts_act = Ts_sol[0]

    Q_in = Qgaz 

    Q_out = Qconds(Ts_act, Tb_actuelle) + Qcondsw(Tb_actuelle)

    derivee = (Q_in - Q_out) / (var.m_Cp)
    
    return derivee


Tb_initial = [300] 

t_span = (0, 172800*2.5) # 2 jours en secondes

t_eval = np.linspace(t_span[0], t_span[1], 500)

# Resolution
# method='RK45' est basically ode45
solution = solve_ivp(fun=dTdt, 
                     t_span=t_span, 
                     y0=Tb_initial, 
                     method='RK45', 
                     t_eval=t_eval)

# Graphique

plt.figure(figsize=(10, 6))
plt.plot(solution.t, solution.y[0], label='Température du corps ($T_b$)', color='red', linewidth=2)

plt.title('Évolution de la température du robot en régime transitoire', fontsize=14)
plt.xlabel('Temps (secondes)', fontsize=12)
plt.ylabel('Kelvin (K)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()

plt.show()