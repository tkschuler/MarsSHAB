import math
import numpy as np
import mars_main
from termcolor import colored
from math import radians
import matplotlib.pyplot as plt
import mars_radiation

import config_mars

'''
m = mars_radiation.MarsRadiation()
h_ang = np.arange(start=5, stop=19, step=.1)

Ls = radians(config.mars_properties['Ls'])
lat = radians(config.mars_properties['lat'])

G_h = []
G_bh = []
G_dh = []
for h in np.nditer(h_ang):
    G_h.append(m.get_global_irradiance(Ls,lat,h))
    G_bh.append(m.get_beam_irradiance(Ls,lat,h))
    G_dh.append(m.get_global_irradiance(Ls,lat,h)-m.get_beam_irradiance(Ls,lat,h))
plt.figure(1)
plt.plot(h_ang,G_h)
plt.plot(h_ang,G_bh)
plt.plot(h_ang,G_dh)
plt.xlabel('Time (hr)')
plt.ylabel('Solar Intensity (W/m^2)')


x = np.linspace(0, 24*3600, 3600)
y = np.linspace(0, 20000, 100)
x_axis, y_axis = np.meshgrid(x, y)
z= np.zeros((100,3600))

direct = []
albedo = []
IR = []

for i in range(0,x.size):
    for j in range(0,y.size):
        h = x[i]
        el = y[j]
        total_radiation = m.get_global_irradiance(Ls,lat,h) + m.get_albedo_flux(Ls,lat,h) + m.get_Mars_IR(el)
        #print h, el, total_radiation
        z[j,i] = total_radiation

for i in range(0,x.size):
    h = x[i]
    direct.append(m.get_global_irradiance(Ls,lat,h))
    albedo.append(m.get_albedo_flux(Ls,lat,h))
    IR.append(m.get_Mars_IR(0))

z_min, z_max = np.amin(z), np.amax(z)

fig, ax = plt.subplots()

x = x/3600
y = y/1000

c = ax.pcolormesh(x, y, z, cmap='magma', vmin=z_min, vmax=z_max)
ax.set_title('Total Available Solar Radiation (W/m^2)')
# set the limits of the plot to the limits of the data
ax.axis([x.min(), x.max(), y.min(), y.max()])
ax.set_xlabel('Time (hr)')
ax.set_ylabel('Elevation (Km)')
fig.colorbar(c, ax=ax)

plt.title('Total Solar Intensity at Martian Surface Level')
#plt.plot(x,direct, label="Direct Solar Flux")
#plt.plot(x,albedo, label="Reflected Albedo Flux")
#plt.plot(x,IR, label="IR Radiation from Surface")
plt.xlabel('Time (hr)')
plt.ylabel('Elevation (km)')
#plt.show()
'''

print("----------------------------------NO MASS FLOW------------------------------------------")

T_s = [210.]
T_i = [210.]
el =[0.]
v= [0.]
vent=0
timespace = np.linspace(0*3600, 24*3600, 24*3600)
alt_dz = 2000.
vel_dz = .1

m = mars_main.Mars_Main()
T_s,T_i,el,v = m.solve_states(T_s,T_i,el,v,timespace,alt_dz,vel_dz)
plt.figure(3,figsize=(10, 8))
plt.plot(timespace/3600,el)
plt.xlabel('Time (hr)')
plt.ylabel('Elevation (m)')

plt.figure(5,figsize=(10, 8))
plt.plot(timespace/(3600),T_s,label="Surface Temperature")
plt.plot(timespace/(3600),T_i,label="Internal Temperature")
plt.xlabel('Time (hr)')
plt.ylabel('Temperature (K)')
plt.title('Solar Balloon Temperature - Mars')
plt.legend(loc='upper right')

''' Comparing Mass Flow Rates'''
'''
print "----------------------------------MASS FLOW RATE 2g/s------------------------------------------"
T_s2 = [210.]
T_i2 = [210.]
el2 =[0.]
v2 = [0.]
mp = 5
vent2 =.006
alt_dz = 1980
vel_dz = -.05

m = mars_main.Mars_Main()
T_s,T_i,el,v = m.solve_states(T_s2,T_i2,el2,v2,mp,vent2,timespace,alt_dz,vel_dz)

plt.figure(6,figsize=(10, 8))
plt.plot(timespace/(3600),T_s,label="Surface Temperature")
plt.plot(timespace/(3600),T_i,label="Internal Temperature")
plt.xlabel('Time (hr)')
plt.ylabel('Temperature (K)')
plt.legend(loc='upper right')

plt.figure(7,figsize=(10, 8))
plt.plot(timespace/(3600),el,label="Surface Temperature")
plt.xlabel('Time (hr)')
plt.ylabel('Elevation (m)')
plt.legend(loc='upper right')
'''
plt.show()

'''
print "----------------------------------MASS FLOW RATE 4g/s------------------------------------------"
T_s2 = [210.]
T_i2 = [210.]
el2 =[0.]
v2= [0.]
mp = 2
vent=.004
timespace = np.linspace(0*3600, 24*3600, 24*3600)
alt_dz = 2040
vel_dz = -.05

T_s2,T_i2,el2,v2 = m.solve_states(T_s2,T_i2,el2,v2,mp,vent,timespace,alt_dz,vel_dz)
print "----------------------------------MASS FLOW RATE 6g/s------------------------------------------"
T_s3 = [210.]
T_i3 = [210.]
el3 =[0.]
v3= [0.]
mp = 2
vent=.006
timespace = np.linspace(0*3600, 24*3600, 24*3600)
alt_dz = 2040
vel_dz = -.05

T_s3,T_i3,el3,v3 = m.solve_states(T_s3,T_i3,el3,v3,mp,vent,timespace,alt_dz,vel_dz)



SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


plt.figure(1,figsize=(10, 8))
plt.plot(timespace/(3600),T_s,label="Surface Temperature")
plt.plot(timespace/(3600),T_i,label="Internal Temperature")
plt.xlabel('Time (hr)')
plt.ylabel('Temperature (K)')
plt.legend(loc='upper right')

plt.figure(2,figsize=(10, 8))
plt.plot(timespace/(3600),el, label="$\dot{m}$ = 2g/s)")
plt.plot(timespace/(3600),el2, label="$\dot{m}$ = 4g/s)")
plt.plot(timespace/(3600),el3, label="$\dot{m}$ = 6g/s)")
plt.xlabel('Time (hr)')
plt.ylabel('Elevation (m)')

plt.axhline(y=1985, color='C0', linestyle='--')
plt.axhline(y=2020, color='C1', linestyle='--')
plt.axhline(y=2040, color='C2', linestyle='--')
plt.axhline(y=2000, color='r', linestyle='-', label = 'Altitude Setpoint')
plt.legend(loc='upper right')

plt.show()
'''
