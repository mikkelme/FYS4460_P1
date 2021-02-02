import numpy as np
import matplotlib.pyplot as plt
# from read_dump import read_dump_ovito
# import subprocess
import plot_set


from b import read_energy



dt_read = [0.02, 0.01, 0.008, 0.002]
energy_files = []
for dt in dt_read:
    dt_str = str(dt).split(".")[-1]
    energy_files.append("etotal" + dt_str + ".txt")

# Plot temperature over time
fig_Temp1 = plt.figure(num = 1, dpi=80, facecolor='w', edgecolor='k')
k = 1; N = 4000
for energy_file in energy_files:
    data = read_energy(energy_file)
    start_idx = len(data)//4
    dt = float("0." + energy_file.strip("etotal.txt"))
    Temp = 2/3*data[:,1]/(N*k)
    Avg_temp = np.mean(Temp[start_idx:])
    std = np.std(Temp[start_idx:])
    plt.plot(data[start_idx:,0]*dt, Temp[start_idx:], label = f"dt = {dt}" + r", $\langle T \rangle$ = " + f"{Avg_temp:.4f}, " + r"$\sigma$ = " + f"{std:.4f}")
plt.xlabel(r"$t/\tau$", fontsize = 14)
plt.ylabel(r"$T/\epsilon k_B^{-1}$ ", fontsize = 14)
plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
plt.legend(loc = "best", fontsize = 13)
fig_Temp1.savefig("../article/figures/temp.pdf", bbox_inches="tight")



#Plot temperature over time depending on system size
energy_files = ["etot_size10.txt", "etot_size12.txt", "etot_size14.txt", "etot_size16.txt"]
system_size = [10, 12, 14, 16]
k = 1
fig_Temp2 = plt.figure(num = 2, dpi=80, facecolor='w', edgecolor='k')
for i in range(len(energy_files)):
    energy_file = energy_files[i]
    data = read_energy(energy_file)
    start_idx = len(data)//4
    size = system_size[i]
    N = 4*size**3
    Temp = 2/3*data[:,1]/(N*k)
    std = np.std(Temp[start_idx:])
    plt.plot(data[start_idx:,0]*dt, Temp[start_idx:], label = f"Size: {size}x{size}x{size}, N = {N}, " + r"$\sigma$ = " + f"{std:.4f}")
plt.xlabel(r"$t/\tau$", fontsize = 14)
plt.ylabel(r"$P/\epsilon \sigma^{-3}$", fontsize = 14)
plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
plt.legend(fontsixe = 13)
fig_Temp2.savefig("../article/figures/temp_size.pdf", bbox_inches="tight")
plt.show()







#
