import numpy as np
import matplotlib.pyplot as plt
from read_dump import read_dump_ovito
import subprocess
import plot_set


def read_energy(energy_file):
    lines = []
    with open(energy_file, "r") as infile:
        info = infile.readline() + infile.readline()
        for line in infile:
            lines.append(np.array(line.split(), dtype=float))
    return np.array(lines)

def run_simulations(runtime, dt):
    script_file = "b.in"
    setting_file = "b_run_settings.in"
    dt_str = str(dt).split(".")[-1]
    with open(setting_file, "w") as outfile:
        outfile.write(f"variable run equal {int(runtime/dt)}\
        \nvariable timestep equal {dt}\
        \nvariable outfile string etotal{dt_str}.txt")
    subprocess.run(["lmp_serial < " +  script_file], shell = True)


if __name__ == "__main__":
    Run = False
    Read = True

    runtime = 20
    dt_run = [0.02, 0.01, 0.008, 0.002]
    dt_read = [0.02, 0.01, 0.008, 0.002]


    if Run:
        for dt in dt_run:
            run_simulations(runtime, dt)
    if Read:
        energy_files = []
        for dt in dt_read:
            dt_str = str(dt).split(".")[-1]
            energy_files.append("etotal" + dt_str + ".txt")

        # Plot total energy over time
        fig_Etot = plt.figure(num = 1, dpi=80, facecolor='w', edgecolor='k')
        for energy_file in energy_files:
            data = read_energy(energy_file)
            dt = float("0." + energy_file.strip("etotal.txt"))
            plt.plot(data[:,0]*dt, data[:,3], label = f"dt = {dt}")

        plt.xlabel(r"$t/\tau$", fontsize = 14)
        plt.ylabel("$U/\epsilon$", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        fig_Etot.savefig("../article/figures/Etot.pdf", bbox_inches="tight")
        plt.legend()
        plt.show()
