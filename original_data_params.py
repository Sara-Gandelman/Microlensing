from data_load import t_cut, F_cut, F_error_cut, I, error
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def chi_square(f_obs, f_exp, f_obs_error):
    chi_sq = sum(((f_obs[i] - f_exp[i]) ** 2) / ((f_obs_error[i]) ** 2) for i in range(len(f_obs)))
    return chi_sq / len(f_obs)


def parabola(x, a, b, c):
    return a * (x ** 2) + b * x + c


# fitting the original data to a parabola without errors, to get an estimation for the parabola parameters
# blabla = np.polyfit(t_cut, F_cut, 2)
blabla = np.polyfit(t_cut, F_cut, 2)

a_init = blabla[0]
b_init = blabla[1]
c_init = blabla[2]

# fitting the parabola with errors
coefficients = curve_fit(parabola, t_cut, F_cut, sigma=F_error_cut, p0=[a_init, b_init, c_init], absolute_sigma=True,
                         xtol=10 ** -8, maxfev=5000)

# coefficients for parabolic fit a * x**2 + b * x + c

a = coefficients[0][0]
b = coefficients[0][1]
c = coefficients[0][2]

f_exp = (t_cut ** 2 * a + t_cut * b + c)

a_error = np.sqrt(abs(coefficients[1][0][0]))
b_error = np.sqrt(abs(coefficients[1][1][1]))
c_error = np.sqrt(abs(coefficients[1][2][2]))

# plot1 = plt.subplot()
# plot2 = plt.subplot()
# plot1.plot(t_cut, f_exp, label='parabola fit')
# plot2.plot(t_cut, F_cut, 'r.', label='observed values')
# plt.legend()
# plt.show()


# T max is the t value for which the flux is maximal, f max is the maximal flux
# f_star is the star's natural flux without the lensing
# miu max is the maximal amount of lensing and u_min is beta(t)/theta einstein

T_max = -b / (2 * a)
f_max = c - (b ** 2) / (4 * a)
f_star = 10 ** (-I[0] / 2.5) * 2550
f_star_error = f_star * np.log(10) * error[0] / 2.5
miu_max = f_max / f_star
u_min = np.sqrt((4 * miu_max ** 2 - 4 - 4 * miu_max * np.sqrt(miu_max ** 2 - 1)) / (2 - 2 * miu_max ** 2))
mu = F_cut / f_star
u = np.sqrt((4 * mu ** 2 - 4 - 4 * mu * np.sqrt(mu ** 2 - 1)) / (2 - 2 * mu ** 2))

T_max_error = np.sqrt((b_error / (2 * a)) ** 2 + (-b * a_error / (2 * a ** 2)) ** 2)
f_max_error = np.sqrt(c_error ** 2 + (2 * b * b_error / (4 * a)) ** 2 + (b ** 2 * a_error / (4 * a ** 2)) ** 2)
miu_max_error = f_max_error / f_star
u_min_error = miu_max_error / (
        (miu_max ** 2 - 1) ** (3 / 2) * np.sqrt(2 * miu_max / (np.sqrt(miu_max ** 2 - 1)) - 2))

print("u_min = ", u_min, "+-", u_min_error)
print("T_max = ", T_max, "+-", T_max_error)
print("f_star = ", f_star, "+-", f_star_error)

t_ranges = np.linspace(561, 621, 1000)
plt.plot(t_cut, F_cut, 'r.', label="original data")
plt.plot(t_ranges, a * t_ranges ** 2 + b * t_ranges + c, label="fitted function")
plt.xlabel("t[HJD]")
plt.ylabel("f[Jy]")
plt.title(" flux from original data vs flux from fitted function as a function of time")
plt.legend()
plt.show()

print(chi_square(F_cut, f_exp, F_error_cut))
