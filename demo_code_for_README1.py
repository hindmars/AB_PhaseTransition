import Module_fphi_AtoB_V01 as AB
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,1)
phi = np.arange(0., 6.5, 0.1, dtype=np.float64)
for p in range(0, 34, 2): ax.plot(phi, AB.fphi_bar(p, 0.5*(10**(-3)), phi))
ax.set_xlabel(r"$\phi/k_{b} \dot T_{c}$");ax.set_ylabel(r"$f(\phi)/(N(0) \times (kb \dot Tc)^{2})$");ax.grid(True)

plt.show()
