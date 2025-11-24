import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
import shutil

methods = ["explicit_euler", "implicit_euler", "improved_euler", "crank_nicolson"]
timesteps = np.linspace(10, 200, 10, dtype=int)
t_end = 4*np.pi
explicit_euler_larger_endtimes = np.linspace(2*np.pi, 20*np.pi, 10)


for method in methods:
    # remove folder if it already exists
    if os.path.exists(method):
        shutil.rmtree(method)

    os.makedirs(method)

    for timestep in timesteps:

        shell_command = f"../build/test_ode {timestep} {t_end} {method}"
        subprocess.run(shell_command, shell = True)
        print(shell_command)

        data = np.loadtxt('output_test_ode.txt', usecols=(0, 1, 2))

        plt.plot(data[:,0], data[:,1], label='position')
        plt.plot(data[:,0], data[:,2], label='velocity')
        plt.xlabel('time')
        plt.ylabel('value')
        plt.title(f'{method}: Mass-Spring System Time Evolution')
        plt.legend()
        plt.grid()
        plt.savefig(f"{method}/{method}_time_plot_timesteps_{timestep}.png")
        plt.close()

        plt.plot(data[:,1], data[:,2], label='phase plot')
        plt.xlabel('position')
        plt.ylabel('velocity')
        plt.title(f'{method}: Mass-Spring System Phase Plot')
        plt.legend()
        plt.grid()
        plt.savefig(f"{method}/{method}_phase_plot_timesteps_{timestep}.png")
        plt.close()


folder = "explicit_euler_larger_endtimes"
if os.path.exists(folder):
    shutil.rmtree(folder)
os.makedirs(folder)
for endtime in explicit_euler_larger_endtimes:
    shell_command = f"../build/test_ode 100 {endtime} explicit_euler"
    subprocess.run(shell_command, shell = True)
    print(shell_command)

    data = np.loadtxt('output_test_ode.txt', usecols=(0, 1, 2))

    plt.plot(data[:,0], data[:,1], label='position')
    plt.plot(data[:,0], data[:,2], label='velocity')
    plt.xlabel('time')
    plt.ylabel('value')
    plt.title("explicit_euler: Mass-Spring System Time Evolution")
    plt.legend()
    plt.grid()
    plt.savefig(f"{folder}/explicit_euler_time_plot_endtime_{endtime}.png")
    plt.close()

    plt.plot(data[:,1], data[:,2], label='phase plot')
    plt.xlabel('position')
    plt.ylabel('velocity')
    plt.title("explicit_euler: Mass-Spring System Phase Plot")
    plt.legend()
    plt.grid()
    plt.savefig(f"{folder}/explicit_euler_phase_plot_endtime_{endtime}.png")
    plt.close()


