import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import mars_radiation
from math import radians

SB = 5.670373E-8

T_amb = 200
T_s = 200

d = 5.79
abs = .6
emiss = .03

Q_sun = abs*500*0.25*math.pi*d*d

Q_loss = SB*2*math.pi*math.pow(d,2)*math.pow(T_s,4)*emiss

Ls = radians(153)
lat = radians(22.3)
d = 20 #m
el = 0
M = 11.3 #kg
cf = 320 #specific heat of aluminim

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dy/dt

def Q_loss(T_s):
    return SB*math.pi*math.pow(d,2)*math.pow(T_s,4)*emiss#*2


def model(T_s,t):
    m = mars_radiation.MarsRadiation()
    Q_sun = m.get_rad_total(lat,Ls,el,t,d)

    dTdt = (Q_sun-Q_loss(T_s))/(cf*M)
    return dTdt

# initial condition
T_s0 = 200

# time points
t = np.linspace(0,100)

# solve ODE
T_s = odeint(model,T_s0,t)

# plot results
plt.plot(t,T_s)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()

d= 50
SF=500
emiss = .03
abs = .8

#T_s = pow((SF*abs)/(emiss*SB)-pow(200,4),.25)
#T_s = pow((abs*500*0.25*math.pi*d*d)/(emiss*SB*2*math.pi*math.pow(d,2))-pow(200,4),.25)
#print T_s
