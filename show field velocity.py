import h5py as h5
import numpy as np
from matplotlib.patches import Circle 
import matplotlib.pyplot as plt
from pylab import *
from PIL import Image
import matplotlib as m

#cdict = {
#  'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
#  'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
#  'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
#}
#
#cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)



name = 'NO.4_0050.h5'
f = h5.File(name)
nx = f['/Nx'][0]
ny = f['/Ny'][0]
fig = plt.figure()
ax = fig.add_subplot(111)

# draw circle
pos = np.array(f['/Pposition'])
pr = np.array(f['/PR'])
pisfree = np.array(f['PIsFree'])
px = pos[0:-2:3]
py = pos[1:-1:3]
for i in range(np.size(px)):
    if(px[i]>=-5 and px[i]<=nx+5 and py[i]>=-5 and py[i]<=ny+5):
        if(pisfree[i]>0):
            cir = Circle(xy=(px[i],py[i]),radius=pr[i],facecolor='k')
        else:
            cir = Circle(xy=(px[i],py[i]),radius=pr[i],facecolor='#363737')

    ax.add_patch(cir)

vel = np.array(f['/Velocity_0'])
vx = vel[0:-2:3]
vy = vel[1:-1:3]
v = np.sqrt(vx**2+vy**2)
v = 10*v.reshape((ny,nx))
norm = m.colors.Normalize(vmin=0, vmax=0.20)
ax0 = ax.pcolor(v, cmap='jet', norm=norm)

plt.axis('equal')
ax.set_xlim(0,nx)
ax.set_ylim(0,ny)

fig.colorbar(ax0, ax=ax)


plt.savefig('/home/foggy/dlvo-master/NO.4_0050.tif' , dpi=300)
plt.show()

