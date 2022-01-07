from part3funcs import grid_search_4D, F_exp, miu_t, u_t
from part2funcs import F, F_error, f_star, t
from data_load import t2, F2, F2_error, f2_star
import matplotlib.pyplot as plt
import numpy as np

# starting ranges for parameters
u_start = 0
u_end = 1
T_start = 560
T_end = 620
f_bl_start = 0
f_bl_end = 1
tau_start = 20
tau_end = 100

params = grid_search_4D([[u_start, u_end], [T_start, T_end], [f_bl_start, f_bl_end], [tau_start, tau_end]], 0.1, 10,
                        0.1, 10, t, F, F_error, f_star)
chi_sq_arr = params[0]
chi_min_idx = params[14]

u_start, u_end = params[15], params[16]
u_amount = params[17]

T_start, T_end = params[18], params[19]
T_amount = params[20]

f_bl_start, f_bl_end = params[21], params[22]
f_bl_amount = params[23]

tau_start, tau_end = params[24], params[25]
tau_amount = params[26]

u_min = params[1]
T_max = params[4]
f_bl = params[7]
tau = params[10]
print(" u_min = ", u_min, " + ", params[3], " - ", params[2])
print(" T_max = ", T_max, " + ", params[6], " - ", params[5])
print(" f_bl = ", f_bl, " + ", params[9], " - ", params[8])
print(" tau = ", tau, " + ", params[12], " - ", params[11])
print("chi square min = ", params[13])

# t_range_1 = np.linspace(-570, 790, 5000)
# f_fit_1 = F_exp(f_star, f_bl, miu_t(u_t(t_range_1, u_min, T_max, tau)))
#
# plt.plot(t, F, 'r.', label="original data")
# plt.plot(t_range_1, f_fit, label="fitted function")
# plt.xlabel("t[HJD]")
# plt.ylabel("f[Jy]")
# plt.title(" flux from original data vs flux from fitted function as a function of time")
# plt.show()

# t_range_2 = np.linspace(-584, 788, 5000)
# f_fit_2 = F_exp(f2_star, f_bl, miu_t(u_t(t_range_2, u_min, T_max, tau)))
# plt.plot(t2, F2, 'r.', label="original data")
# plt.plot(t_range_2, f_fit_2, label="fitted function")
# plt.xlabel("t[HJD]")
# plt.ylabel("f[Jy]")
# plt.title(" flux from original data vs flux from fitted function as a function of time")
# plt.show()

# t_range_3 = np.linspace(-576, 788, 5000)
# f_fit_3 = F_exp(f_star, f_bl, miu_t(u_t(t_range_3, u_min, T_max, tau)))
# plt.plot(t2, F2, 'r.', label="original data")
# plt.plot(t_range_3, f_fit_3, label="fitted function")
# plt.xlabel("t[HJD]")
# plt.ylabel("f[Jy]")
# plt.title(" flux from original data vs flux from fitted function as a function of time")
# plt.show()

arr_u_T = chi_sq_arr[:, :, chi_min_idx[2], chi_min_idx[3]]

arr_u_f_bl = chi_sq_arr[:, chi_min_idx[1], :, chi_min_idx[3]]

arr_u_tau = chi_sq_arr[:, chi_min_idx[1], chi_min_idx[2], :]

arr_T_f_bl = chi_sq_arr[chi_min_idx[0], :, :, chi_min_idx[3]]

arr_T_tau = chi_sq_arr[chi_min_idx[0], :, chi_min_idx[2], :]

arr_f_bl_tau = chi_sq_arr[chi_min_idx[0], chi_min_idx[1], :, :]

u_range = np.linspace(u_start, u_end, u_amount)
T_range = np.linspace(T_start, T_end, T_amount)
f_bl_range = np.linspace(f_bl_start, f_bl_end, f_bl_amount)
tau_range = np.linspace(tau_start, tau_end, tau_amount)

X, Y = np.meshgrid(u_range, T_range)
fig1, ax = plt.subplots(1)
arr_u_T = arr_u_T.T
b = ax.contourf(X, Y, arr_u_T)
ax.plot(u_min, T_max, 'r.', label="chi_square minimum")
plt.legend()
plt.xlabel("u_min")
plt.ylabel("T_0[HJD]")
plt.xlim([u_start, u_end])
plt.ylim([T_start, T_end])
plt.title("u_min and T_0 chi square contour plot")
plt.clabel(b, inline=1, fontsize=15, colors='r')
plt.show()

X, Y = np.meshgrid(u_range, f_bl_range)
fig2, ax = plt.subplots(1)
arr_u_f_bl = arr_u_f_bl.T
b = ax.contourf(X, Y, arr_u_f_bl)
ax.plot(u_min, f_bl, 'r.', label="chi_square minimum")
plt.legend()
plt.xlabel("u_min")
plt.ylabel("f_bl")
plt.xlim([u_start, u_end])
plt.ylim([f_bl_start, f_bl_end])
plt.title("u_min and f_bl chi square contour plot")
plt.clabel(b, inline=1, fontsize=15, colors='r')
plt.show()

X, Y = np.meshgrid(u_range, tau_range)
fig3, ax = plt.subplots(1)
arr_u_tau = arr_u_tau.T
b = ax.contourf(X, Y, arr_u_tau)
ax.plot(u_min, tau, 'r.', label="chi_square minimum")
plt.legend()
plt.xlabel("u_min")
plt.ylabel("tau[HJD]")
plt.xlim([u_start, u_end])
plt.ylim([tau_start, tau_end])
plt.title("u_min and tau chi square contour plot")
plt.clabel(b, inline=1, fontsize=15, colors='r')
plt.show()

X, Y = np.meshgrid(T_range, f_bl_range)
fig4, ax = plt.subplots(1)
arr_T_f_bl = arr_T_f_bl.T
b = ax.contourf(X, Y, arr_T_f_bl)
ax.plot(T_max, f_bl, 'r.', label="chi_square minimum")
plt.legend()
plt.xlabel("T_0[HJD]")
plt.ylabel("f_bl")
plt.ylim([f_bl_start, f_bl_end])
plt.xlim([T_start, T_end])
plt.title("T_0 and f_bl chi square contour plot")
plt.clabel(b, inline=1, fontsize=15, colors='r')
plt.show()

X, Y = np.meshgrid(T_range, tau_range)
fig5, ax = plt.subplots(1)
arr_T_tau = arr_T_tau.T
b = ax.contourf(X, Y, arr_T_tau)
ax.plot(T_max, tau, 'r.', label="chi_square minimum")
plt.legend()
plt.xlabel("T_0[HJD]")
plt.ylabel("tau[HJD]")
plt.ylim([tau_start, tau_end])
plt.xlim([T_start, T_end])
plt.title("T_0 and tau chi square contour plot")
plt.clabel(b, inline=1, fontsize=15, colors='r')
plt.show()

X, Y = np.meshgrid(f_bl_range, tau_range)
fig6, ax = plt.subplots(1)
arr_f_bl_tau = arr_f_bl_tau.T
b = ax.contourf(X, Y, arr_f_bl_tau)
ax.plot(f_bl, tau, 'r.', label="chi_square minimum")
plt.legend()
plt.ylabel("tau[HJD]")
plt.xlabel("f_bl")
plt.xlim([f_bl_start, f_bl_end])
plt.ylim([tau_start, tau_end])
plt.title("f_bl and tau chi square contour plot")
plt.clabel(b, inline=1, fontsize=15, colors='r')
plt.show()
