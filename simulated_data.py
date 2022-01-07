from data_load import t_cut, F_cut, F_error_cut, I
import numpy as np
from scipy.optimize import curve_fit
from original_data_params import parabola, f_star, T_max
import matplotlib.pyplot as plt
from random import gauss

T_max_arr = []
f_max_arr = []
u_min_arr = []
tau_arr = []

T_max_error_arr = []
f_max_error_arr = []
u_min_error_arr = []
tau_error_arr = []

for i in range(10 ** 4):
    t_simul = np.random.choice(t_cut, len(t_cut))
    F_simul = []
    for k in range(len(t_simul)):
        F_simul.append(float(F_cut[t_cut == t_simul[k]]))
    F_simul = np.array(F_simul)
    sigma = F_simul / 10
    F_error_simul = gauss(sigma, sigma / 10)

    init_coefficients = np.polyfit(t_simul, F_simul, 2)

    # coefficients for parabolic fit a * x**2 + b * x + c
    a_init_sim = init_coefficients[0]
    b_init_sim = init_coefficients[1]
    c_init_sim = init_coefficients[2]

    coefficients = curve_fit(parabola, t_simul, F_simul, sigma=F_error_simul, p0=[a_init_sim, b_init_sim, c_init_sim])

    a_sim = coefficients[0][0]
    b_sim = coefficients[0][1]
    c_sim = coefficients[0][2]

    T_max_sim = -b_sim / (2 * a_sim)
    f_max_sim = c_sim - (b_sim ** 2) / (4 * a_sim)
    miu_max_sim = f_max_sim / f_star
    if miu_max_sim < 1:
        continue
    u_min_sim_squared = (4 * miu_max_sim ** 2 - 4 - 4 * miu_max_sim * np.sqrt(miu_max_sim ** 2 - 1)) / (
            2 - 2 * miu_max_sim ** 2)
    u_min_sim = np.sqrt(u_min_sim_squared)

    T_max_arr.append(T_max_sim)
    f_max_arr.append(f_max_sim)
    u_min_arr.append(u_min_sim)

