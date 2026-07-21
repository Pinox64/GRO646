import numpy as np

# Géométrie et propriétés du robot
D = 0.10  # Diamètre externe du robot [m]
Db = 0.05  # Diamètre interne (corps) [m]
A_cont = 1e-4  # Surface de contact [m^2]
A_s = 4 * np.pi * (D/2)**2 - A_cont  # Surface extérieure [m^2]
x_ins = 0.025  # Épaisseur isolant [m]
k_in = 0.005  # Conductivité isolant [W/(m.K)]
rho_Cp_b = 3e6  # Capacité calorifique volumique du corps [J/(m^3.K)]
V_b = (4/3) * np.pi * (Db/2)**3  # Volume du corps [m^3]
m_Cp = rho_Cp_b * V_b  # Capacité thermique massique totale du corps [J/K]

# Chaleur générée
Q_fc = 1.0  # Puissance de combustion [W]
eta_tot = 0.02  # Rendement
Q_gen = (1 - eta_tot) * Q_fc  # Chaleur réelle dissipée par le corps [W]

# Rayonnement et surfaces
eps_s = 0.01  # Émissivité sans VEC
eps_vec = 0.09  # Émissivité VEC
A_vec = 0.75 * A_s  # Surface VEC [m^2]
A_rad = A_s - A_vec  # Surface sans VEC [m^2]
sigma = 5.67e-8  # Constante de Stefan-Boltzmann [W/(m^2.K^4)]

# Commutateur thermique (Switch)
k_sw = 59.0  # Conductivité du switch [W/(m.K)]
x_sw = 0.03  # Longueur du switch [m]
A_sw = 1.5e-6  # Surface de contact du switch [m^2]

# Environnement (Mars)
g = 3.71  # Gravité [m/s^2]
T_sun = 5780.0  # Température du Soleil [K]
r_s = 0.696e9  # Rayon du Soleil [m]