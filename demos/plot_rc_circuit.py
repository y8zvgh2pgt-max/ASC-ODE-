import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../build/output_rc.txt")
t = data[:, 0]
Uc = data[:, 1]
U0 = data[:, 2]

plt.figure(figsize=(10, 6))
plt.plot(t, U0, label="Quelle U0(t) (Input)", linestyle="--", alpha=0.6)
plt.plot(t, Uc, label="Kondensator Uc(t) (Output)", linewidth=2)

plt.xlabel("Zeit t [s]")
plt.ylabel("Spannung [V]")
plt.title("RC-Kreis Simulation (Crank-Nicolson)")
plt.legend()
plt.grid()
plt.savefig("rc_circuit_plot.png")
