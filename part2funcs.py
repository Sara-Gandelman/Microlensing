import sys
import time
from data_load import t, I, error
from original_data_params import f_star
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare

tau = 60.799
tau_sigma = 0.041
F = 10 ** (-I / 2.5) * 2550
F_error = F * np.log(10) * 0.003 / 2.5
mu = F / f_star
mu_error = F_error / f_star


def chi_square(f_obs, f_exp, f_obs_error):
    chi_sq = sum(((f_obs[i] - f_exp[i]) ** 2) / ((f_obs_error[i]) ** 2) for i in range(len(f_obs)))
    return chi_sq / len(f_obs)


def u_t(t, u_min, t_0, tau):
    return np.sqrt(u_min ** 2 + ((t - t_0) / tau) ** 2)


def miu_t(u):
    return (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))


def F_exp(f_star, miu):
    return f_star * miu


def find_error(array, min_array):
    i = min_array[0]
    j = min_array[1]

    # finding the distance of 1 sigma for 2nd parameter while keeping the 1st parameter fixed
    right_movement_1 = 0
    left_movement_1 = 0
    right_movement_2 = 0
    left_movement_2 = 0
    for k in range(j, len(array[0])):
        if array[i][k] - array[i][j] >= 2.3:
            right_movement_2 = k - j
            break

    for k in range(j, -1, -1):
        if array[i][k] - array[i][j] >= 2.3:
            left_movement_2 = j - k
            break

    # finding the distance of 1 sigma for 1st parameter while keeping the 2nd parameter fixed

    for k in range(i, len(array)):
        if array[k][j] - array[i][j] >= 2.3:
            right_movement_1 = k - i
            break

    for k in range(i, -1, -1):
        if array[k][j] - array[i][j] >= 2.3:
            left_movement_1 = i - k
            break

    return left_movement_1, left_movement_2, right_movement_1, right_movement_2


# i values are u_min values and j values are t_v values
def grid_chi_min(parameters_ranges, step_u, step_T, time, f, f_error):
    i_start = parameters_ranges[0][0]
    i_end = parameters_ranges[0][1]
    j_start = parameters_ranges[1][0]
    j_end = parameters_ranges[1][1]
    i_sigma_right = 0
    i_sigma_left = 0
    j_sigma_right = 0
    j_sigma_left = 0
    i_min = 0
    j_min = 0
    while step_u > 10 ** -3:
        chi_square_arr = np.zeros((int((i_end - i_start) / step_u), int((j_end - j_start) / step_T)))
        chi_square_min = sys.float_info.max
        chi_min_idx = []
        k = 0
        for i in np.linspace(i_start, i_end, int((i_end - i_start) / step_u)):
            m = 0
            for j in np.linspace(j_start, j_end, int((j_end - j_start) / step_T)):
                f_exp = F_exp(f_star, miu_t(u_t(time, i, j, tau)))
                chi_sq = chi_square(f, f_exp, f_error)
                # if chi_sq < 9.71:
                #     plt.plot(2, 4)
                #     plt.show()
                chi_square_arr[k, m] = chi_sq
                if chi_sq < chi_square_min:
                    chi_square_min = chi_sq
                    i_min = i
                    j_min = j
                    chi_min_idx = [k, m]
                m += 1
            k += 1

        if step_u * 0.1 <= 10 ** -3:
            break
        i_start = i_min - step_u if i_min - step_u >= i_start else i_start
        i_end = i_min + step_u if i_min + step_u < i_end else i_end
        j_start = j_min - step_T if j_min - step_T > j_start else j_start
        j_end = j_min + step_T if j_min + step_T < j_end else j_end
        step_u = step_u * 0.1
        step_T = step_T * 0.1

        # finding the error of the best parameters using the chi squared grid

        left_movement_1, left_movement_2, right_movement_1, right_movement_2 = find_error(chi_square_arr,
                                                                                          chi_min_idx)
        # now we know the amount of movement left and right for both axis, and so we can calculate the 1 sigma error as
        # the amount of movement for each side * the step size
        step_size_i = (i_end - i_start) / m
        step_size_j = (j_end - j_start) / k

        i_sigma_right = step_size_i * right_movement_1 if right_movement_1 != 0 else i_sigma_right
        i_sigma_left = step_size_i * left_movement_1 if left_movement_1 != 0 else i_sigma_left
        j_sigma_right = step_size_j * right_movement_2 if right_movement_2 != 0 else j_sigma_right
        j_sigma_left = step_size_j * left_movement_2 if left_movement_2 != 0 else j_sigma_left

        # print(chi_square_arr)

    return chi_square_arr, i_min, i_sigma_left, i_sigma_right, j_min, j_sigma_left, j_sigma_right, chi_square_min, chi_min_idx, i_start, i_end, j_start, j_end, step_u, step_T


# t1 = time.time_ns()
# params = grid_chi_min([[0, 1], [560, 620]], 0.1, 10, t, F, F_error)
# u_min = params[1]
# T_max = params[4]
# min_idx = params[8]
# t2 = time.time_ns()
# chi_square_arr = params[0]
# print(" u_min = ", params[1], " + ", params[3], " - ", params[2])
# print(" t_0 = ", params[4], " + ", params[6], " - ", params[5])
# print("chi square min = ", params[7])
# print("time taken: ", (t2 - t1) / 10 ** 9)

# we need this!!!!!!!!!!!!! uncomment at the end!!!!!
# t_range = np.linspace(-570, 790, 5000)
# plt.plot(t, F, 'r.', label="original data")
# plt.plot(t_range, F_exp(f_star, miu_t(u_t(t_range, params[1], params[4], tau))), label="fitted function")
# plt.xlabel("t[HJD]")
# plt.ylabel("f[Jy]")
# plt.title(" flux from original data vs flux from fitted function as a function of time")
# plt.legend()
# plt.show()

# u_min_start, u_min_end, t_0_start, t_0_end, step_u, step_T = params[9], params[10], params[11], params[12], params[13], \
#                                                              params[14]
#
# u_min_range = np.linspace(u_min_start, u_min_end, 20, endpoint=True)
# t_0_range = np.linspace(t_0_start, t_0_end, 20, endpoint=True)
# print(u_min_range)
# print(t_0_range)
# #
# chi_square_arr = chi_square_arr.T
# X, Y = np.meshgrid(u_min_range, t_0_range)
# fig1, ax = plt.subplots(1)
# b = ax.contourf(X, Y, chi_square_arr)
# # ax.imshow(chi_square_arr)
# ax.plot(u_min, T_max, '.', label="chi_square minimum")
# # ax.plot(min_idx[1], min_idx[0], '.', label="chi_square minimum")
# plt.legend()
# plt.title("chi square contour map of parameters u_min and T_max")
# plt.xlabel("u_min")
# plt.ylabel("T_max[HJD]")
# plt.ylim([t_0_start, t_0_end])
# plt.xlim([u_min_start, u_min_end])
# plt.clabel(b, inline=1, fontsize=15)
# plt.show()

