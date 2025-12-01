import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("output_pendulum.txt")
t     = data[:, 0]
theta = data[:, 1]
omega = data[:, 2]

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(t, theta)
plt.title("pendulum angle")
plt.xlabel("t")
plt.ylabel(r"$\theta(t)$")
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(t, omega)
plt.title("pendulum angular velocity")
plt.xlabel("t")
plt.ylabel(r"$\dot{\theta}(t)$")
plt.grid()

plt.tight_layout()
plt.savefig("pendulum_plot.png")
