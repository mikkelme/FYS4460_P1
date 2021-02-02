import numpy as np
import matplotlib.pyplot as plt
import subprocess
import plot_set
import statsmodels.api as sm

def run_simulations(temp, density ):
    script_file = "e.in"
    setting_file = "e_run_settings.in"
    T_str = str(temp).replace(".","")
    d_str = str(density).replace(".","")
    with open(setting_file, "w") as outfile:
        outfile.write(f"variable temp equal {temp}\
        \nvariable density equal {density}\
        \nvariable outfile string e_datafile_T{T_str}_d{d_str}.txt")
    subprocess.run(["lmp_serial < " +  script_file], shell = True)


if __name__ == "__main__":
    Run = False
    Read = True

    temp_run = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    density_run = [0.01, 0.02, 0.03, 0.04, 0.05]
    dt = 0.002
    size = 10

    from b import read_energy

    if Run:
        for temp in temp_run:
            for density in density_run:
                run_simulations(temp, density)


    if Read:
        data_files = []
        for temp in temp_run:
            for density in density_run:
                T_str = str(temp).replace(".","")
                d_str = str(density).replace(".","")
                data_files.append(f"e_datafile_T{T_str}_d{d_str}.txt")


        avg_temp = []
        density = []
        avg_pressure = []
        for data_file in data_files:
            data = read_energy(data_file)
            start_idx = len(data)//2
            N = 4*size**3
            Temp = 2/3*data[start_idx:,1]/(N)
            avg_temp.append(np.mean(Temp))
            avg_pressure.append(np.mean(data[start_idx:,4]))
            density.append(float("0." + data_file.split("d")[-1].strip(".txt")[-2:]))

        avg_temp = np.array(avg_temp)
        density = np.array(density)
        avg_pressure = np.array(avg_pressure)
        Txd = avg_temp*density


        figd = plt.figure(num=1, dpi=80, facecolor='w', edgecolor='k')
        temp_space = np.linspace(avg_temp[0], avg_temp[-1], int(1e3))
        for i in range(5):
            idx = np.argwhere(density == density[i])
            plot = plt.plot(avg_temp[idx], avg_pressure[idx], "o", label = f"d = {density[i]} " + r"$\rho/\sigma^{-3}$")
            x,y = avg_temp[idx], avg_pressure[idx]
            x = sm.add_constant(x)
            model = sm.OLS(y, x)
            res = model.fit()
            b, a = res.params
            b_err, a_err = res.bse
            plt.plot(temp_space, temp_space*a + b, color = plot[0].get_color(), linestyle ="--", label = f"Lin. fit: a = {a:.6f}" + r" $\pm$ " + f"{a_err:.1g}")

        plt.xlabel(r"$T/\epsilon k_B^{-1}$ ", fontsize = 14)
        plt.ylabel(r"$P/\epsilon \sigma^{-3}$", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
        plt.savefig("../article/figures/temp_press.pdf", bbox_inches="tight")


        figTxd = plt.figure(num=2, dpi=80, facecolor='w', edgecolor='k')
        x,y = Txd, avg_pressure
        x = sm.add_constant(x)
        model = sm.OLS(y, x)
        res = model.fit()
        b, a = res.params
        b_err, a_err = res.bse
        Txd_space = np.linspace(np.min(Txd), np.max(Txd), int(1e3))
        plt.plot(Txd, avg_pressure, "o", label = r"$T \cdot \rho$ datapoints ")
        plt.plot(Txd_space, Txd_space*a + b, linestyle ="--", label = f"Linear fit:\na = {a:.2f}" + r" $\pm$ " + f"{a_err:.1g}\nb = {b:.4f}" + r" $\pm$ " + f"{b_err:.1g}")

        plt.xlabel(r"$T\rho/\epsilon k_B^{-1} \sigma^{-3}$", fontsize = 14)
        plt.ylabel(r"$P/\epsilon \sigma^{-3}$", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
        plt.savefig("../article/figures/temp_rho_press.pdf", bbox_inches="tight")


        plt.show()





        #plt.show()
