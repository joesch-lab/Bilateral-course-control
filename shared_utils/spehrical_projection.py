import numpy as np
# import custom_projection_example
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = np.sqrt(XsqPlusYsq + z**2)               # r
    elev = np.atan2(z, np.sqrt(XsqPlusYsq))      # theta
    az = np.atan2(y,x)                           # phi
    return r, elev, az

def sph2cart(az, elev, r):
    z = r * np.sin(elev)
    rcoselev = r * np.cos(elev)
    x = rcoselev * np.cos(az)
    y = rcoselev * np.sin(az)
    return x, y, z

def get_xyz(radius):
    x, y = np.mgrid[-radius:radius, -radius:radius]
    center = [0, 0, 0]
    r = np.sqrt((center[0] - x) ** 2 + (center[1] - y) ** 2)
    z = np.sqrt(radius ** 2 - r ** 2)
    return x, y, z

def sinusoidal_projection(az, elev):
    x = az * np.cos(elev)
    y = elev
    return x, y

def mollweide_projection(az, elev):
    map = Basemap(projection='mol',lat_0=0,lon_0=0,resolution='l')
    x, y = map(az, elev)
    return x, y

projection = "mollweide"
fig1, ax1 = plt.subplots(subplot_kw={'projection': projection})
ax1.set_longitude_grid_ends(90)
# fig2, ax2 = plt.subplots()
# img = plt.imread("pinwheel.png")
radius = 50

total_az = 0
total_elev = 0
total_az, total_elev = np.mgrid[0:2*np.pi:360, np.pi:np.pi:180]
##front and back
ax1.plot([-np.pi/2, -np.pi/2], [-np.pi/2, np.pi/2], color='k')
ax1.plot([np.pi/2, np.pi/2], [-np.pi/2, np.pi/2], color='k')
##left eye
ax1.fill_betweenx(y=[-np.pi/2, np.pi/2], x1=-np.pi + (20/180)*np.pi, x2=(10/180)*np.pi, color='grey', alpha=0.5, label='Total', lw=0)
##right eye
ax1.fill_betweenx(y=[-np.pi/2, np.pi/2], x1=-(10/180)*np.pi, x2=np.pi - (20/180)*np.pi, color='grey', alpha=0.5, lw=0)
##ephys
ax1.fill_betweenx(y=[(25/90) * -np.pi/2, (65/90) * np.pi/2], x1=(70/90)*-np.pi/2, x2=(70/90)*np.pi/2, color='g', alpha=0.2, label='Electrophysiology', lw=0)
##behaviour
ax1.fill_betweenx(y=[(10/90) * np.pi/2, np.pi/2], x1=-np.pi, x2=np.pi, color='r', alpha=0.2, label='Behaviour', lw=0)
ax1.legend(loc='upper right')
ax1.grid(visible=True, linewidth=0.5, color='grey', alpha=0.5)
ax1.spines.clear()
plt.show()
fig1.savefig(projection + 'projection.svg', dpi=600, transparent=True)
fig1.savefig(projection + 'projection.pdf', dpi=600, transparent=True)
fig1.savefig(projection + 'projection.png', dpi=300)

# epsilon = 0.0001
# elev_range = [0 + (10/90)*np.pi/2, np.pi/2] #-pi to pi, 0 is at the horizon
# az_range = [0, 2*np.pi + epsilon] #0 to 2*pi
# az, elev = np.mgrid[az_range[0]:az_range[1]:(az_range[1]-az_range[0])/128, elev_range[0]:elev_range[1]:(elev_range[1]-elev_range[0])/128]
# x, y, z = sph2cart(az, elev, 50)
# # x ,y, z = get_xyz(radius)
# print(x.shape)
# print(np.amax(az))
# facecolors = np.sin(az*6) + 1
# # facecolors = facecolors.astype('uint8')
# # print(facecolors.shape, np.amax(facecolors))
# facecolors = np.dstack((facecolors, facecolors, facecolors, facecolors))
# ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=img, alpha=0.2)
#
# ## draw a circle that denotes a fly
# az, elav = np.mgrid[-np.pi:np.pi:2*np.pi/50, -np.pi:np.pi:2*np.pi/50]
# x, y, z = sph2cart(az, elav, 5)
# # x ,y, z = get_xyz(radius)
# ax.plot_surface(x, y, z, rstride=2, cstride=2, cmap='Reds')
#
# ax.set_xlim3d(-50, 50)
# ax.set_ylim3d(-50, 50)
# ax.set_zlim3d(-50, 50)
# ax.view_init(elev=45, azim=45)
# plt.show()
# fig.savefig('example.png')