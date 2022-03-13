


import numpy as np
import matplotlib.pyplot as plt

# box = dict(facecolor='yellow', pad=5, alpha=0.2)

# fig = plt.subplots(2, 2)
# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
# fig.subplots_adjust(left=0.2, wspace=0.6)

# c1 = ax1.scatter(rad, r)

rad = np.array([[0.5, 1.2],[0.9, 3.0]]) * np.pi
r = np.array([[0.5, 0.1],[0.9, 0.3]])

fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2, 2, subplot_kw=dict(projection='polar'))
ax1.scatter(rad, r)
ax2.scatter(rad, r)
ax3.scatter(rad, r)
ax4.scatter(rad, r)

plt.show()

'''
ax1.plot(2000 * np.random.rand(10))
ax1.set_title('ylabels not aligned')
ax1.set_ylabel('misaligned 1', bbox=box)
ax1.set_ylim(0, 2000)

ax3.set_ylabel('misaligned 2', bbox=box)
ax3.plot(np.random.rand(10))

xlabel = -0.3  # axes coords

ax2.set_title('ylabels aligned')
ax2.plot(2000 * np.random.rand(10))
ax2.set_ylabel('aligned 1', bbox=box)
ax2.yaxis.set_label_coords(xlabel, 0.5)
ax2.set_ylim(0, 2000)

ax4.plot(np.random.rand(10))
ax4.set_ylabel('aligned 2', bbox=box)
ax4.yaxis.set_label_coords(xlabel, 0.5)
'''


