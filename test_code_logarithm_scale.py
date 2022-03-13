import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import bivariate_normal
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages


delta = 0.5

x = np.arange(-3.0, 4.001, delta)
y = np.arange(-4.0, 3.001, delta)
X, Y = np.meshgrid(x, y)
Z = bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)

fig  = plt.figure()
ax   = fig.add_subplot(1,1,1)
#axim = ax.imshow(Z, norm = LogNorm())
#axim   = ax.contourf(X,Y,Z,levels=[1e-3,1e-2,1e-1,1e0],cmap=plt.cm.jet,norm = LogNorm())
axim    = ax.contourf(X,Y,Z,20,cmap=plt.cm.jet,norm = LogNorm())
cb   = fig.colorbar(axim)

pp = PdfPages('fig.pdf')
pp.savefig()
pp.close()


plt.show()
