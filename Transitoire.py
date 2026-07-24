import Variables as var
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

temp_choice = 0
wind_choice = 0
VEC_choice = True
Comm_choice = True
h_nat = 0

setups = [
    {"label": "Cas 1: T_inf=50K, vent=0", "temp_choice": 0, "wind_choice": 0},
    {"label": "Cas 2: T_inf=50K, vent=2", "temp_choice": 0, "wind_choice": 2},
    {"label": "Cas 3: T_inf=250K, vent=0", "temp_choice": 4, "wind_choice": 0},
    {"label": "Cas 4: T_inf=250K, vent=2", "temp_choice": 4, "wind_choice": 2},
]

# Constant
Qgaz = 0.98
r_out = var.D/2
r_in = var.Db/2

# Fonction de Ts
Qrad = lambda Ts: var.eps_s * var.sigma * (var.A_s-var.A_vec) * (Ts**4 - var.T[temp_choice]**4)
Qrad_vec = lambda Ts: var.eps_vec * var.sigma * (var.A_vec) * (Ts**4 - var.T[temp_choice]**4)
Qconv = lambda Ts, h_nat: h_nat*var.A_s * (Ts - var.T[temp_choice])

# Fonction de Tb
Qcondsw = lambda Tb: var.k_sw*var.A_sw*(Tb - var.T[temp_choice])/var.x_sw
Qcondg = lambda Tb: var.k_in*var.A_cont*(Tb - var.T[temp_choice])/var.x_ins

# Fonction de Ts et Tb
Qconds = lambda Ts, Tb: (4 * np.pi * var.k_in * r_in * r_out / (r_out - r_in))*(Tb - Ts)

def bilan_surface(Ts_guess, Tb_actuelle):
    Ts = Ts_guess[0] # fsolve passe une liste
    
    chaleur_venant_du_corps = Qconds(Ts, Tb_actuelle)
    
    tf = (Ts + var.T[temp_choice])/2

    Beta = 1/tf

    Ra_D = (var.g * Beta * abs(Ts - var.T[temp_choice]) * (var.D**3) * var.Pr[temp_choice]) / (var.nu[temp_choice]**2)

    if wind_choice == 0:
        Nu_sphere_Nat = 2 + (0.589*(Ra_D**(1/4)))/(1 + (0.469/var.Pr[temp_choice])**(9/16))**(4/9)
        h = (Nu_sphere_Nat * var.k[temp_choice])/var.D

    elif wind_choice == 2:
        #put the wind equation here
        Re=wind_choice*var.D/var.nu[temp_choice]
        u_inf=var.rho[temp_choice]*wind_choice*var.D/Re

        Close_index = np.argmin([np.abs(x - Ts) for x in var.T])
        Re_s=wind_choice*var.D/var.nu[Close_index]
        u_s= var.rho[Close_index]*wind_choice*var.D/Re_s

        Nuforce=2+(0.4*np.sqrt(Re)+0.06*(Re**(2/3)))*(var.Pr[temp_choice]**0.4)*((u_inf/u_s)**0.25)
        h=Nuforce*var.k[temp_choice]/var.D

    # Bilan = 0
    residu = chaleur_venant_du_corps + Qsun - Qrad(Ts) - (Qrad_vec(Ts) if VEC_choice else 0) - Qconv(Ts, h)
    return residu

def dTdt(t, Tb_actuelle):
    Tb_val = Tb_actuelle[0] # fsolve passe une liste
    Ts_sol = fsolve(bilan_surface, x0=[300], args=(Tb_val,))
    Ts_act = Ts_sol[0]

    Q_in = Qgaz 

    Q_out = Qconds(Ts_act, Tb_val) + (Qcondsw(Tb_val) if Comm_choice else 0) + Qcondg(Tb_val)

    derivee = (Q_in - Q_out) / (var.m_Cp)
    
    return [derivee]

Tb_initial = [300] 

t_span = (0, 172800*1.5) # 2 jours en secondes

t_eval = np.linspace(t_span[0], t_span[1], 500)

def solve_transitoire(Tb_initial, t_span, t_eval):
    # Resolution
    # method='RK45' est basically ode45
    solution = solve_ivp(fun=dTdt, 
                        t_span=t_span, 
                        y0=Tb_initial, 
                        method='RK45', 
                        t_eval=t_eval)
    return solution

for config in setups:
    temp_choice = config["temp_choice"]
    wind_choice = config["wind_choice"]
    # Temp changed so Qsun changed too
    Qsun = var.sigma * var.eps_s * (var.A_s / 2) * (var.r_s/(var.r_max if temp_choice == 0 else var.r_min))**2 * var.T_sun**4

    plt.figure(figsize=(10, 6))

    VEC_choice = False
    Comm_choice = False
    solution = solve_transitoire(Tb_initial, t_span, t_eval)
    # Retrouver Tb
    Ts_calcules = []
    for Tb_val in solution.y[0]:
        Ts_sol = fsolve(bilan_surface, x0=[300], args=(Tb_val,))
        Ts_calcules.append(Ts_sol[0])
    plt.plot(solution.t, solution.y[0], label="Tb : VEC - OFF, Comm - OFF", color='orange', linewidth=2)
    plt.plot(solution.t, Ts_calcules, label="Ts : VEC - OFF, Comm - OFF", color='yellow', linewidth=2)

    VEC_choice = True
    Comm_choice = False
    solution = solve_transitoire(Tb_initial, t_span, t_eval)
    # Retrouver Tb
    Ts_calcules = []
    for Tb_val in solution.y[0]:
        Ts_sol = fsolve(bilan_surface, x0=[300], args=(Tb_val,))
        Ts_calcules.append(Ts_sol[0])
    plt.plot(solution.t, solution.y[0], label="Tb : VEC - ON, Comm - OFF", color='blue', linewidth=2)
    plt.plot(solution.t, Ts_calcules, label="Ts : VEC - ON, Comm - OFF", color='cyan', linewidth=2)

    VEC_choice = False
    Comm_choice = True
    solution = solve_transitoire(Tb_initial, t_span, t_eval)
    # Retrouver Tb
    Ts_calcules = []
    for Tb_val in solution.y[0]:
        Ts_sol = fsolve(bilan_surface, x0=[300], args=(Tb_val,))
        Ts_calcules.append(Ts_sol[0])
    plt.plot(solution.t, solution.y[0], label="Tb : VEC - OFF, Comm - ON", color='green', linewidth=2)
    plt.plot(solution.t, Ts_calcules, label="Ts : VEC - OFF, Comm - ON", color='lightgreen', linewidth=2)

    VEC_choice = True
    Comm_choice = True
    solution = solve_transitoire(Tb_initial, t_span, t_eval)
    # Retrouver Tb
    Ts_calcules = []
    for Tb_val in solution.y[0]:
        Ts_sol = fsolve(bilan_surface, x0=[300], args=(Tb_val,))
        Ts_calcules.append(Ts_sol[0])
    plt.plot(solution.t, solution.y[0], label="Tb : VEC - ON, Comm - ON", color='red', linewidth=2)
    plt.plot(solution.t, Ts_calcules, label="Ts : VEC - ON, Comm - ON", color='magenta', linewidth=2)

    plt.axhline(y=350, color='red', linestyle='--', linewidth=1.5, label='Limite Max (350 K)')
    plt.axhline(y=200, color='blue', linestyle='--', linewidth=1.5, label='Limite Min (200 K)')
    plt.title(config["label"], fontsize=14)
    plt.xlabel('Temps (secondes)', fontsize=12)
    plt.ylabel('Kelvin (K)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()

plt.show()