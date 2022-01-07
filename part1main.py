import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from original_data_params import T_max, f_max, u_min, T_max_error, f_max_error, u_min_error
from simulated_data import T_max_arr, f_max_arr, u_min_arr, tau_arr, T_max_error_arr, f_max_error_arr, u_min_error_arr, \
    tau_error_arr

plt.hist(T_max_error_arr, np.linspace(0, 250, 100))
plt.show()

n_T, bins_T, patches_T = plt.hist(T_max_arr, np.linspace(T_max - 2, T_max + 2, 300))
plt.show()

n_f, bins_f, patches_f = plt.hist(f_max_arr, np.linspace(0.0066, 0.0068, 300))
plt.show()

n_u, bins_u, patches_u = plt.hist(u_min_arr, np.linspace(0.67, 0.7, 300))
plt.show()

n_tau, bins_tau, patches_tau = plt.hist(tau_arr, np.linspace(0, 1000, 200))
plt.xlim([0, 1450])
plt.ylim([0, 122])
plt.show()


def gaussian(x, a, sigma, miu):
    return a * np.e ** ((-(x - miu) ** 2) / (2 * sigma ** 2))


T_max_fit_parameters = curve_fit(gaussian, bins_T[:len(bins_T) - 1], n_T, p0=[810, 3, 590])

T_a = T_max_fit_parameters[0][0]
T_sigma = T_max_fit_parameters[0][1]
T_miu = T_max_fit_parameters[0][2]

T_a_error = np.sqrt(abs(T_max_fit_parameters[1][0][0]))
T_sigma_error = np.sqrt(abs(T_max_fit_parameters[1][1][1]))
T_miu_error = np.sqrt(abs(T_max_fit_parameters[1][2][2]))

plt.hist(T_max_arr, np.linspace(T_max - 2, T_max + 2, 300))
plt.plot(bins_T, gaussian(bins_T, T_a, T_sigma, T_miu))
plt.title(" T_max gaussian fit vs T_max histogram")
plt.xlabel("t [HJD]")
plt.ylabel("count")
plt.xlim([T_miu - 5, T_miu + 5])
plt.show()

f_max_fit_parameters = curve_fit(gaussian, bins_f[:len(bins_f) - 1], n_f, p0=[500, 3 * 10 ** -4, 0.0064])

f_a = f_max_fit_parameters[0][0]
f_sigma = f_max_fit_parameters[0][1]
f_miu = f_max_fit_parameters[0][2]
#
plt.hist(f_max_arr, np.linspace(0.0066, 0.0068, 300))
plt.plot(bins_f, gaussian(bins_f, f_a, f_sigma, f_miu))
plt.title(" f_max gaussian fit vs f_max histogram")
plt.xlabel("f [Jy]")
plt.ylabel("count")
plt.xlim([0.0065, 0.00678])
plt.show()

f_a_error = np.sqrt(abs(f_max_fit_parameters[1][0][0]))
f_sigma_error = np.sqrt(abs(f_max_fit_parameters[1][1][1]))
f_miu_error = np.sqrt(abs(f_max_fit_parameters[1][2][2]))

u_min_fit_parameters = curve_fit(gaussian, bins_u[:len(bins_u) - 1], n_u, p0=[250, 0.02, 0.69])

u_a = u_min_fit_parameters[0][0]
u_sigma = u_min_fit_parameters[0][1]
u_miu = u_min_fit_parameters[0][2]

plt.hist(u_min_arr, np.linspace(0.67, 0.7, 300))
plt.plot(bins_u, gaussian(bins_u, u_a, u_sigma, u_miu))
plt.title(" u_min gaussian fit vs u_min histogram")
plt.xlabel("u_min")
plt.ylabel("count")
plt.show()

u_a_error = np.sqrt(abs(u_min_fit_parameters[1][0][0]))
u_sigma_error = np.sqrt(abs(u_min_fit_parameters[1][1][1]))
u_miu_error = np.sqrt(abs(u_min_fit_parameters[1][2][2]))

tau_fit_parameters = curve_fit(gaussian, bins_tau[:len(bins_tau) - 1], n_tau, p0=[122, 50, 200])

tau_a = tau_fit_parameters[0][0]
tau_sigma = tau_fit_parameters[0][1]
tau_miu = tau_fit_parameters[0][2]

plt.hist(tau_arr, np.linspace(0, 1000, 200))
plt.plot(bins_tau, gaussian(bins_tau, tau_a, tau_sigma, tau_miu))
plt.title(" tau gaussian fit vs u_min histogram")
plt.xlabel("t [HJD]")
plt.ylabel("count")
plt.xlim([0, 1450])
plt.show()

T_a_error = np.sqrt(abs(T_max_fit_parameters[1][0][0]))
T_sigma_error = np.sqrt(abs(T_max_fit_parameters[1][1][1]))
T_miu_error = np.sqrt(abs(T_max_fit_parameters[1][2][2]))


print(" T_0 from the original data came out to be : ", T_max, " +- ", T_max_error)
print(" T_0 from the simulated data came out to be : ", T_miu, " +- ", T_sigma)

print("\n f_max from the original data came out to be : ", f_max, " +- ", f_max_error)
print(" f_max from the simulated data came out to be : ", f_miu, " +- ", f_sigma)

print("\n u_min from the original data came out to be : ", u_min, " +- ", u_min_error)
print(" u_min from the simulated data came out to be : ", u_miu, " +- ", u_sigma)

