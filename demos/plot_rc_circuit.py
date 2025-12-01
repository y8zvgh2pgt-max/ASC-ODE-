import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import subprocess


methods = ["explicit_euler", "implicit_euler", "improved_euler", "crank_nicolson"]
timesteps = np.linspace(10, 200, 10, dtype=int)
main_dir_name = "rc_circuit"


for method in methods:
    folder = f"{main_dir_name}/{method}"

    # remove folder if it already exists
    if os.path.exists(f"{folder}"):
        shutil.rmtree(f"{folder}")

    os.makedirs(f"{folder}")

    for timestep in timesteps:

        shell_command = f"../build/demo_rc_circuit {timestep} {method}"
        subprocess.run(shell_command, shell = True)
        print(shell_command)

        data = np.loadtxt('output_rc.txt', usecols=(0, 1, 2))

        plt.plot(data[:,0], data[:,2], label="source U0(t) (Input)", linestyle="--", alpha=0.6)
        plt.plot(data[:,0], data[:,1], label="capacitor Uc(t) (Output)", linewidth=2)
        plt.xlabel("time t [s]")
        plt.ylabel("voltage [V]")
        plt.title(f"{method}: rc-circuit simulation for {timestep} timesteps.")
        plt.legend()
        plt.grid()
        plt.savefig(f"{folder}/{method}_rc_circuit_plot_timesteps_{timestep}.png")
        plt.close()
