import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(0, 100)
ax.set_ylim(0, 30)
ax.set_zlim(-10, 0)
ax.bar3d(0, 0, 0, 20, 10, -1, alpha=0.5,
         color='lightcoral')
ax.bar3d(0, 0, -2, 20, 10, -2, color='r', alpha=0.5)
ax.bar3d(0, 20, 0, 20, 10, -1, color='lightcoral', alpha=0.5)
ax.bar3d(0, 20, -2, 20, 10, -2, color='r', alpha=0.5)
ax.bar3d(30, 5, 0, 20, 20, -1, color='lightcoral', alpha=0.5)
ax.bar3d(30, 5, -2, 20, 20, -2, color='r', alpha=0.5)
ax.bar3d(30, 5, -5, 20, 20, -5, color='darkred', alpha=0.5)
ax.bar3d(0, 0, -5, 20, 10, -5, color='darkred', alpha=0.5)
ax.bar3d(20, 5, -5, 10, 5, -5, color='darkred', alpha=0.5)
plt.show()
