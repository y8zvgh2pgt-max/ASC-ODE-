import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("output_legendre.txt", skiprows=1)
x = data[:, 0]

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
for n in range(6): # order 0 to 5
    col_idx = 1 + 2*n 
    plt.plot(x, data[:, col_idx], label=f'$P_{n}(x)$')

plt.title("legendre polynomials")
plt.xlabel("x")
plt.grid()
plt.legend()

plt.subplot(1, 2, 2)
for n in range(1, 6):
    col_idx = 1 + 2*n + 1
    plt.plot(x, data[:, col_idx], label=f"$P'_{n}(x)$")

plt.title("derivatives (via AutoDiff)")
plt.xlabel("x")
plt.grid()
plt.legend()

plt.tight_layout()
plt.savefig("legendre_plot.png")