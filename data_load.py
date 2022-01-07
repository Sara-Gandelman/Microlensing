import numpy as np
import matplotlib.pyplot as plt

f = open("2019-BLG-0035.txt", "r")
t = []
I = []
error = []

while True:
    line = f.readline().split()

    if line:
        I.append(float(line[1]))
        t.append(float(line[0]))
        error.append(float(line[2]))
    if not line:
        break

error = np.array(error)
t = np.array(t)
t -= 2458000
I = np.array(I)
# plt.plot(t, I,'b.')
# plt.show()
I_cut = I[I < 14.1]
t_cut = t[I < 14.1]
# plt.plot(t_cut, I_cut, 'b.')
# plt.show()
error_cut = error[I < 14.1]
# F_cut = I_cut
# F_error_cut = error_cut
F_cut = 10 ** (-I_cut / 2.5) * 2550
F_error_cut = F_cut * np.log(10) / 2.5 * error_cut

# making sure everything is fine
# plot = plt.subplot()
# plt.errorbar(t_cut, F_cut, yerr=F_error_cut, fmt=".")
# plt.show()

f = open("f_bl!=1.txt", "r")
t2 = []
I2 = []
error2 = []

while True:
    line = f.readline().split()

    if line:
        I2.append(float(line[1]))
        t2.append(float(line[0]))
        error2.append(float(line[2]))
    if not line:
        break

error2 = np.array(error2)
t2 = np.array(t2)
t2 -= 2458000
I2 = np.array(I2)
F2 = 10 ** (-I2 / 2.5) * 2550
F2_error = F2 * np.log(10) / 2.5 * error2
f2_star = 10 ** (-I2[0] / 2.5) * 2550
# plt.plot(t2, I2, 'b.')
# plt.show()
