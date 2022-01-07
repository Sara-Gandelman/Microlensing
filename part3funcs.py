import sys
from part2funcs import miu_t, u_t, chi_square, t, find_error, F, F_error
from original_data_params import f_star
from data_load import F2, t2, F2_error, f2_star
import numpy as np
import matplotlib.pyplot as plt


def find_error_4D(array, min_array):
    u_idx = min_array[0]
    T_idx = min_array[1]
    f_bl_idx = min_array[2]
    tau_idx = min_array[3]

    # finding the distance of 1 sigma for 2nd parameter while keeping the 1st parameter fixed
    right_movement_u = 0
    left_movement_u = 0
    right_movement_T = 0
    left_movement_T = 0
    right_movement_f_bl = 0
    left_movement_f_bl = 0
    right_movement_tau = 0
    left_movement_tau = 0

    # u sigmas
    for k in range(u_idx + 1, len(array[:, T_idx, f_bl_idx, tau_idx])):
        if array[k, T_idx, f_bl_idx, tau_idx] - array[u_idx, T_idx, f_bl_idx, tau_idx] >= 4.72:
            right_movement_u = k - u_idx
            break

    for k in range(u_idx, -1, -1):
        if array[k, T_idx, f_bl_idx, tau_idx] - array[u_idx, T_idx, f_bl_idx, tau_idx] >= 4.72:
            left_movement_u = u_idx - k
            break

    # T sigmas

    for k in range(T_idx + 1, len(array[u_idx, :, f_bl_idx, tau_idx])):
        if array[u_idx, k, f_bl_idx, tau_idx] - array[u_idx, T_idx, f_bl_idx, tau_idx] >= 4.72:
            right_movement_T = k - T_idx
            break

    for k in range(T_idx, -1, -1):
        if array[u_idx, k, f_bl_idx, tau_idx] - array[u_idx, T_idx, f_bl_idx, tau_idx] >= 4.72:
            left_movement_T = T_idx - k
            break

    # f_bl sigmas
    for k in range(len(array[u_idx, T_idx, :, tau_idx])):
        if array[u_idx, T_idx, k, tau_idx] - array[u_idx, T_idx, f_bl_idx, tau_idx] >= 4.72:
            right_movement_f_bl = f_bl_idx - k
            break

    for k in range(f_bl_idx, -1, -1):
        if array[u_idx, T_idx, k, tau_idx] - array[u_idx, T_idx, f_bl_idx, tau_idx] >= 4.72:
            left_movement_f_bl = k - f_bl_idx
            break

    # tau sigmas

    for k in range(len(array[u_idx, T_idx, f_bl_idx, :])):
        if array[u_idx, T_idx, f_bl_idx, k] - array[u_idx, T_idx, f_bl_idx, tau_idx] >= 4.72:
            right_movement_tau = k - tau_idx
            break

    for k in range(tau_idx, -1, -1):
        if array[u_idx, T_idx, f_bl_idx, k] - array[u_idx, T_idx, f_bl_idx, tau_idx] >= 4.72:
            left_movement_tau = k - tau_idx
            break

    return left_movement_u, right_movement_u, left_movement_T, right_movement_T, left_movement_f_bl, right_movement_f_bl, left_movement_tau, right_movement_tau


def F_exp(f_star, f_bl, miu):
    return f_star * f_bl * miu + f_star * (1 - f_bl)


def grid_search_4D(parameters_ranges, step_u, step_T, step_f_bl, step_tau, t, f, f_error, f_star):
    u_start = parameters_ranges[0][0]
    u_end = parameters_ranges[0][1]
    T_start = parameters_ranges[1][0]
    T_end = parameters_ranges[1][1]
    f_bl_start = parameters_ranges[2][0]
    f_bl_end = parameters_ranges[2][1]
    tau_start = parameters_ranges[3][0]
    tau_end = parameters_ranges[3][1]

    u_sigma_right = 0
    u_sigma_left = 0
    T_sigma_right = 0
    T_sigma_left = 0
    f_bl_sigma_right = 0
    f_bl_sigma_left = 0
    tau_sigma_right = 0
    tau_sigma_left = 0

    u_min_min = 0
    T_0_min = 0
    f_bl_min = 0
    tau_min = 0

    while step_u > 10 ** -3:
        chi_square_arr = np.zeros((int((u_end - u_start) / step_u), int((T_end - T_start) / step_T),
                                   int((f_bl_end - f_bl_start) / step_f_bl), int((tau_end - tau_start) / step_tau)))
        chi_square_min = sys.float_info.max
        chi_min_idx = []
        k = 0
        for u in np.linspace(u_start, u_end, int((u_end - u_start) / step_u)):
            print(k)
            m = 0
            for T in np.linspace(T_start, T_end, int((T_end - T_start) / step_T)):
                a = 0
                for f_bl in np.linspace(f_bl_start, f_bl_end, int((f_bl_end - f_bl_start) / step_f_bl)):
                    b = 0
                    for tau in np.linspace(tau_start, tau_end, int((tau_end - tau_start) / step_tau)):
                        f_exp = F_exp(f_star, f_bl, miu_t(u_t(t, u, T, tau)))
                        chi_sq = chi_square(f, f_exp, f_error)
                        chi_square_arr[k, m, a, b] = chi_sq
                        if chi_sq < chi_square_min:
                            chi_square_min = chi_sq
                            u_min_min = u
                            T_0_min = T
                            f_bl_min = f_bl
                            tau_min = tau
                            chi_min_idx = [k, m, a, b]
                        b += 1
                    a += 1
                m += 1
            k += 1

        if step_u * 0.1 <= 10 ** -3:
            break
        u_start = u_min_min - step_u if u_min_min - step_u >= u_start else u_start
        u_end = u_min_min + step_u if u_min_min + step_u < u_end else u_end
        T_start = T_0_min - step_T if T_0_min - step_T > T_start else T_start
        T_end = T_0_min + step_T if T_0_min + step_T < T_end else T_end
        f_bl_start = f_bl_min - step_f_bl if f_bl_min - step_f_bl >= f_bl_start else f_bl_start
        f_bl_end = f_bl_min + step_f_bl if f_bl_min + step_f_bl < f_bl_end else f_bl_end
        tau_start = tau_min - step_tau if tau_min - step_tau > tau_start else tau_start
        tau_end = tau_min + step_tau if tau_min + step_tau < tau_end else tau_end

        step_u = step_u * 0.1
        step_T = step_T * 0.1
        step_f_bl = step_f_bl * 0.1
        step_tau = step_tau * 0.1

        # finding the error of the best parameters using the chi squared grid

        left_movement_u, right_movement_u, left_movement_T, right_movement_T, left_movement_f_bl, right_movement_f_bl, left_movement_tau, right_movement_tau = find_error_4D(
            chi_square_arr,
            chi_min_idx)

        # now we know the amount of movement left and right for all parameters,
        # and so we can calculate the 1 sigma error as the amount of movement for each side * the step size

        step_size_u = (u_end - u_start) / k
        step_size_T = (T_end - T_start) / m
        step_size_f_bl = (f_bl_end - f_bl_start) / a
        step_size_tau = (tau_end - tau_start) / b

        u_sigma_right = step_size_u * right_movement_u if right_movement_u != 0 else u_sigma_right
        u_sigma_left = step_size_u * left_movement_u if left_movement_u != 0 else u_sigma_left
        T_sigma_right = step_size_T * right_movement_T if right_movement_T != 0 else T_sigma_right
        T_sigma_left = step_size_T * left_movement_T if left_movement_T != 0 else T_sigma_left
        f_bl_sigma_right = step_size_f_bl * right_movement_f_bl if right_movement_f_bl != 0 else f_bl_sigma_right
        f_bl_sigma_left = step_size_f_bl * left_movement_f_bl if left_movement_f_bl != 0 else f_bl_sigma_left
        tau_sigma_right = step_size_tau * right_movement_tau if right_movement_tau != 0 else tau_sigma_right
        tau_sigma_left = step_size_tau * left_movement_tau if left_movement_tau != 0 else tau_sigma_left

    return chi_square_arr, u_min_min, u_sigma_left, u_sigma_right, T_0_min, T_sigma_left, T_sigma_right, \
           f_bl_min, f_bl_sigma_left, f_bl_sigma_right, tau_min, tau_sigma_left, tau_sigma_right, chi_square_min, chi_min_idx, \
           u_start, u_end, k, T_start, T_end, m, f_bl_start, f_bl_end, a, tau_start, tau_end, b
