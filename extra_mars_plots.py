import math
import numpy as np
import mars_main
from termcolor import colored
from math import radians
import matplotlib.pyplot as plt


print "----------------------------------NO MASS FLOW------------------------------------------"

T_s = [210.]
T_i = [210.]
el =[0.]
v= [0.]
mp =10
vent=0
timespace = np.linspace(0*3600, 48*3600, 48*3600)
alt_dz = 0
vel_dz = 0

m = mars_main.Mars_Main()
T_s,T_i,el,v = m.solve_states(T_s,T_i,el,v,mp,vent,timespace,alt_dz,vel_dz)
plt.figure(3,figsize=(10, 8))
plt.plot(timespace/3600,el)
plt.xlabel('Time (hr)')
plt.ylabel('Elevation (m)')

plt.figure(5,figsize=(10, 8))
plt.plot(timespace/3600,T_s,label="Surface Temperature")
plt.plot(timespace/3600,T_i,label="Internal Temperature")
plt.xlabel('Time (hr)')
plt.ylabel('Temperature (K)')
plt.legend(loc='upper right')

''' Comparing Mass Flow Rates'''

print "----------------------------------MASS FLOW RATE 2g/s------------------------------------------"
T_s = [210.]
T_i = [210.]
el =[0.]
v= [0.]
mp = 2
vent=.015
timespace = np.linspace(0*3600, 24*3600, 24*3600)
alt_dz = 5000
vel_dz = -.01

m = mars_main.Mars_Main()
T_s,T_i,el,v = m.solve_states(T_s,T_i,el,v,mp,vent,timespace,alt_dz,vel_dz)

print "----------------------------------MASS FLOW RATE 4g/s------------------------------------------"
T_s2 = [210.]
T_i2 = [210.]
el2 =[0.]
v2= [0.]
mp =5
vent=.01
timespace = np.linspace(0*3600, 24*3600, 24*3600)
alt_dz = 5010
vel_dz = -.01

T_s2,T_i2,el2,v2 = m.solve_states(T_s2,T_i2,el2,v2,mp,vent,timespace,alt_dz,vel_dz)
print "----------------------------------MASS FLOW RATE 6g/s------------------------------------------"
T_s3 = [210.]
T_i3 = [210.]
el3 =[0.]
v3= [0.]
mp =10
vent=.006
timespace = np.linspace(0*3600, 24*3600, 24*3600)
alt_dz = 5000
vel_dz = -.05

T_s3,T_i3,el3,v3 = m.solve_states(T_s3,T_i3,el3,v3,mp,vent,timespace,alt_dz,vel_dz)

"""PLOTTING"""

'''
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
'''

plt.figure(1,figsize=(10, 8))
plt.plot(timespace/3600,T_s,label="Surface Temperature")
plt.plot(timespace/3600,T_i,label="Internal Temperature")
plt.xlabel('Time (hr)')
plt.ylabel('Temperature (K)')
plt.legend(loc='upper right')

plt.figure(2,figsize=(10, 8))
plt.plot(timespace/3600,el, label="$\dot{m}$ = 2g/s)")
plt.plot(timespace/3600,el2, label="$\dot{m}$ = 4g/s)")
plt.plot(timespace/3600,el3, label="$\dot{m}$ = 6g/s)")
plt.xlabel('Time (hr)')
plt.ylabel('Elevation (m)')

plt.axhline(y=1985, color='C0', linestyle='--')
plt.axhline(y=2020, color='C1', linestyle='--')
plt.axhline(y=2040, color='C2', linestyle='--')
plt.axhline(y=2000, color='r', linestyle='-', label = 'Altitude Setpoint')
plt.legend(loc='upper right')

plt.show()
