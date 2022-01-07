from random import gauss
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from part1main import gaussian
from part2funcs import F, t, tau, grid_chi_min, F_error, F_exp, miu_t, u_t, chi_square
from data_load import I

T_max_arr = []
u_min_arr = []
for i in range(10 ** 4):
    print(i)
    t_simul = t
    F_simul = F + F * np.random.randn(len(F)) / 10
    # t_simul = np.random.choice(t, len(t))
    # F_simul = []
    # for k in range(len(t_simul)):
    #     F_simul.append(float(F[t == t_simul[k]]))
    F_simul = np.array(F_simul)
    sigma_sim = F_simul / 10
    F_error_simul = gauss(sigma_sim, sigma_sim / 10)
    params = grid_chi_min([[0, 1], [560, 620]], 0.1, 10, t_simul, F_simul, F_error_simul)
    u_min_sim, T_0_sim = params[1], params[4]
    u_min_arr.append(u_min_sim)
    T_max_arr.append(T_0_sim)

n_T_max, bins_T_max, patches_T_max = plt.hist(T_max_arr, np.linspace(587, 592, 20))
plt.show()

n_u_min, bins_u_min, patches_u_min = plt.hist(u_min_arr, np.linspace(0.675, 0.715, 20))
plt.show()

T_0_fit_parameters = curve_fit(gaussian, bins_T_max[:len(bins_T_max) - 1], n_T_max, p0=[3000, 3, 590], maxfev=5000)

T_0_a = T_0_fit_parameters[0][0]
T_0_sigma = T_0_fit_parameters[0][1]
T_0_miu = T_0_fit_parameters[0][2]

u_min_fit_parameters = curve_fit(gaussian, bins_u_min[:len(bins_u_min) - 1], n_u_min, p0=[3000, 0.01, 0.69],
                                 maxfev=5000)

u_min_a = u_min_fit_parameters[0][0]
u_min_sigma = u_min_fit_parameters[0][1]
u_min_miu = u_min_fit_parameters[0][2]

print("\nT_0 = ", T_0_miu, "+-", T_0_sigma)
print("u_min = ", u_min_miu, "+-", u_min_sigma)
print("T_0_a = ", T_0_a)
print("u_min_a = ", u_min_a)

plt.hist(T_max_arr, np.linspace(587, 592, 20))
plt.plot(bins_T_max, gaussian(bins_T_max, T_0_a, T_0_sigma, T_0_miu))
plt.title(" T_0 gaussian fit vs T_0 histogram using simulations")
plt.xlabel("t [HJD]")
plt.ylabel("count")
plt.show()

plt.hist(u_min_arr, np.linspace(0.675, 0.715, 20))
plt.plot(bins_u_min, gaussian(bins_u_min, u_min_a, u_min_sigma, u_min_miu))
plt.title(" u_min gaussian fit vs u_min histogram using simulations")
plt.xlabel("u_min")
plt.ylabel("count")
plt.show()

print("T_0_arr = ", T_max_arr)
print("u_min_arr = ", u_min_arr)

f_T_exp = gaussian(bins_T_max, T_0_a, T_0_sigma, T_0_miu)
f_u_exp = gaussian(bins_u_min, u_min_a, u_min_sigma, u_min_miu)

chi_square_u_min_fit = chi_square(n_u_min, f_u_exp,
                                  np.ones((len(t))) * (bins_u_min[1] - bins_u_min[0]) / np.sqrt(12))
chi_square_T_0_fit = chi_square(n_T_max, f_T_exp, np.ones((len(t))) * (bins_T_max[1] - bins_T_max[0]) / np.sqrt(12))

print("chi square for u min fit is : ", chi_square_u_min_fit)
print("chi square for T0 fit is : ", chi_square_T_0_fit)
