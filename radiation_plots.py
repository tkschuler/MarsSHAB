import radiation
import matplotlib.pyplot as plt
import math
import fluids
import numpy as np

# Test Variables

doy = 306 #temporary day of year
lat = math.radians(35.106766) # rad
h_ang = 0
el = 0 #elevation (m)

'''PLOTS -------------------------------------------------------------------------'''


r = radiation.Radiation(doy,lat,h_ang,el)

h_ang = np.arange(start=-90, stop=90, step=.1)

zen2 = []
h_ang2 = []
I_direct_solar = []
I_diffuse_sky=[]
I_reflected_surface = []
for h in np.nditer(h_ang):
    zen2.append(math.pi/2 - r.get_zenith(lat,math.radians(h)))
    h_ang2.append(1./15.*h)

    z = r.get_zenith(lat,math.radians(h))
    I_direct_solar.append(r.get_direct_SI(z,el))
    I_diffuse_sky.append(r.get_diffuse_SI(z,el))
    I_reflected_surface.append(r.get_reflected_SI(z,el))

plt.figure(1)
plt.plot(h_ang2,zen2)
plt.xlabel('Time (hr)')
plt.ylabel('Solar Elevation (deg)')

plt.figure(2)
plt.plot(h_ang2,I_direct_solar)
plt.plot(h_ang2,I_diffuse_sky)
plt.plot(h_ang2,I_reflected_surface)
plt.xlabel('Time (hr)')
plt.ylabel('Solar Intensity (W/m^2)')


I_IR_EARTH = []
I_IR_SKY = []
el = np.arange(start= 0, stop=40000, step=10)
for e in np.nditer(el):
    I_IR_EARTH.append(r.get_earth_IR(e))
    I_IR_SKY.append(r.get_sky_IR(e))

plt.figure(3)
plt.plot(I_IR_EARTH,el)
plt.plot(I_IR_SKY,el)
plt.xlabel('Solar Intensity (W/m^2)')
plt.ylabel('Elevation (m)')

''' Radiation Colormesh '''

x = np.linspace(-6, 6, 100)
y = np.linspace(0, 40000, 100)
x_axis, y_axis = np.meshgrid(x, y)
z= np.zeros((100,100))

for i in range(0,x.size):
    for j in range(0,y.size):
        h = x[i]
        el = y[j]
        #zen = r.get_zenith(lat,math.radians(h*15))
        #I = r.get_direct_SI(zen,el) + r.get_diffuse_SI(zen,el) + r.get_reflected_SI(zen,el)
        #IR = r.get_earth_IR(el)+r.get_sky_IR(el)
        #total_radiation = I+IR

        total_radiation = r.get_rad_total(lat,el,math.radians(h*15),5.79)/100

        print h, el, total_radiation
        z[j,i] = total_radiation

z_min, z_max = np.amin(z), np.amax(z)
#z_min, z_max = 0, 600

fig, ax = plt.subplots()

c = ax.pcolormesh(x, y, z, cmap='magma', vmin=z_min, vmax=z_max)
ax.set_title('Radiation Intensity')
# set the limits of the plot to the limits of the data
ax.axis([x.min(), x.max(), y.min(), y.max()])
ax.set_xlabel('Solar Intensity (W/m^2)')
ax.set_ylabel('Elevation (m)')
fig.colorbar(c, ax=ax)

print "\nTEST"

rad = r.get_rad_total(lat,40000,0,5.79)
print "total radiation:" , rad

plt.show()
