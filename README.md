# AB_PhaseTransition

Repository of modules and codes developed for studying the AB phase transition in superfluid Helium-3, developed in collaboration with GitHub user timohyva. It started as a fork of https://github.com/timohyva/AB_PhaseTransition/.


## Descriptions 

The main modules of the repository are


## Languages (current)
* python3

## License
[GPL] (https://www.gnu.org/licenses/old-licenses/gpl-2.0-faq.en.html)

## Demo

```python
import Module_SC_Beta_V04 as SCC
import Module_AB_wall_V00 as WB
import he3_tools as h
import Module_fphi_AtoB_V01 as AB

import numpy as np
import matplotlib.pyplot as plt

fig1, ax1 = plt.subplots(1,1)
phi = np.arange(0., 6.5, 0.1, dtype=np.float64)
for p in range(0, 34, 2): ax1.plot(phi, AB.fphi_bar(p, 0.5*(10**(-3)), phi))
ax1.set_xlabel(r"$\phi/k_{b} \dot T_{c}$");ax1.set_ylabel(r"$f(\phi)/(N(0) \times (kb \dot Tc)^{2})$");ax1.grid(True)

fig2, ax2 = plt.subplots(1,1)
p1 = np.arange(0, 34.2, 0.2); p2 = np.arange(21.2, 34.2, 0.2)
ax2.plot(h.Tc_mK_expt(p1), p1,color="blue")
ax2.plot(h.TAB_mK_expt(p2), p2,color="purple")
ax2.set_xlabel(r"$T/mK$");ax2.set_ylabel(r"$p/bar$")
ax2.grid(True)

plt.show()
```
