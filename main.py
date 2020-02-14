import math
import fluids
import numpy as np
import radiation
import sphere_balloon
from termcolor import colored
from scipy.integrate import odeint
import matplotlib.pyplot as plt


Rsp_air = 287.058
Cp_air0 = 1003.8
d = 5.79
vol = math.pi*4/3*pow((d/2),3)
cs_area = math.pi*d*d/4.0
surfArea = math.pi*d*d
vm_coeff = .5

massEnv = surfArea*939*7.87E-6#rhoEnv*envThickness
massPayload = 0
cp = 2250.0 #(J/(kg K)) specific heat of envelope material
k = massEnv*cp #thermal mass coefficient


doy = 306 #temporary day of year
lat = math.radians(35.106766) # rad

class State:
    def __init__(self):
        self.v= 0.0
        self.el = 132.6 #m
        atm = fluids.atmosphere.ATMOSPHERE_1976(self.el)
        self.Ti = atm.T #287.288117979
        self.Ts = atm.T #287.288117979 #atm temperature at 132.6 m

class StateDerrivative:
    def __init__(self):
        self.dv = 0.0
        self.dz = 0.0
        self.dTi = 0.0
        self.dTs = 0.0

'''
def get_acceleration(state):
    atm = fluids.atmosphere.ATMOSPHERE_1976(state.el)
    rho_int = atm.P/(Rsp_air*state.Ti)
    rho_atm = atm.rho
    #print "VOL:", vol
    F_b = (rho_atm - rho_int)*vol*atm.g
    Cd = .8
    F_d =  Cd*(0.5*rho_atm*math.fabs(state.v)*state.v)*cs_area#Sphere_Balloon.get_Cd(state.v, state.el)*(0.5*rho_atm*fabs(state.v)*state.v)*cs_area
    vm = (massEnv + massPayload) + rho_atm*vol + vm_coeff*rho_atm*vol
    #print "!!!", massEnv+ massPayload, rho_atm*vol, vm_coeff*rho_atm*vol
    accel = (F_b  - F_d - (massEnv + massPayload)*atm.g)/vm
    #print "T_s", state.Ts, "rho_int:", rho_int, "rho_atm:", rho_atm, "F_b:", F_b, "F_d", F_d,"vm:",vm, "g:", atm.g
    #print accel
    return accel
'''

def get_acceleration(v,el,T_s,T_i):
    atm = fluids.atmosphere.ATMOSPHERE_1976(el)
    rho_int = atm.P/(Rsp_air*T_i)
    rho_atm = atm.rho
    #print "VOL:", vol
    F_b = (rho_atm - rho_int)*vol*atm.g
    Cd = .8
    F_d =  Cd*(0.5*rho_atm*math.fabs(v)*v)*cs_area#Sphere_Balloon.get_Cd(state.v, state.el)*(0.5*rho_atm*fabs(state.v)*state.v)*cs_area
    vm = (massEnv + massPayload) + rho_atm*vol + vm_coeff*rho_atm*vol
    #print "!!!", massEnv+ massPayload, rho_atm*vol, vm_coeff*rho_atm*vol
    accel = (F_b  - F_d - (massEnv + massPayload)*atm.g)/vm
    #print "T_s", T_s, "rho_int:", rho_int, "rho_atm:", rho_atm, "F_b:", F_b, "F_d", F_d,"vm:",vm, "g:", atm.g
    #print accel
    return accel


#rate of change of Surface Temperature
def get_dTs(state,h):
    bal = sphere_balloon.Sphere_Balloon(5.79,0.8)
    rad = radiation.Radiation(doy,lat,h,state.el)
    q_rad  = rad.get_rad_total(lat, state.el, h,5.79)
    q_surf = bal.get_sum_q_surf(q_rad, state.Ts, state.el, state.v)
    q_int  = bal.get_sum_q_int(state.Ts, state.Ti, state.el)
    dT_sdt = (q_surf-q_int)/k
    #print "q_rad: ", q_rad, "q_surf: ", q_surf, "q_int: ", q_int
    #print "dT_sdt:", dT_sdt, "\n"
    return dT_sdt

def get_dTi(state):
    bal = sphere_balloon.Sphere_Balloon(5.79,0.8)
    atm = fluids.atmosphere.ATMOSPHERE_1976(state.el)
    q_int  = bal.get_sum_q_int(state.Ts, state.Ti, state.el)
    tm_air = atm.rho*vol*Cp_air0
    return q_int/tm_air


#----------------------------------------------------------------------------------------------------------


#am i doing this right

def eval(state, dt, h, der):
    state.el = state.el + der.dz*dt
    state.v  = state.v  + der.dv*dt
    if state.el < 132.6 and state.v < 0:
        #state.v  = 0
        state.el = 132.6
    state.Ti = state.Ti + der.dTi*dt
    state.Ts = state.Ts + der.dTs*dt


    out = StateDerrivative()
    out.dz = state.v
    out.dv = get_acceleration(state)
    out.dTs = get_dTs(state,h + dt*((math.pi/6.)/(3600.))) #dt is converted to hour angle
    out.dTi = get_dTi(state)

    for t1, w1 in zip(t, wsol):
            print t1, w1[0], w1[1], w1[2], w1[3], w1[4]

    return odeint(vectorfield, w0, t,atol=abserr, rtol=relerr)

def integrate(state, h, dt):
    zero = StateDerrivative()

    #These are all state derrivatives
    a = eval(state, 0.0, h, zero)
    b = eval(state, 0.5*dt, h, a)
    c = eval(state, 0.5*dt, h, b)
    d = eval(state, dt, h, c)

    del_dt = (1./6.)*(a.dz + 2*b.dz + 2*c.dz + d.dz);
    dv_dt  = (1./6.)*(a.dv  + 2*b.dv  + 2*c.dv  + d.dv);
    dTi_dt = (1./6.)*(a.dTi + 2*b.dTi + 2*c.dTi + d.dTi);
    dTs_dt = (1./6.)*(a.dTs + 2*b.dTs + 2*c.dTs + d.dTs);

    state.el += del_dt*dt;
    state.v += dv_dt*dt;

    if state.el < 132.6 and state.v < 0: #balloon does not descend below el = min_el
        #state.v  = 0
        state.el = 132.6

    state.Ti += dTi_dt*dt
    state.Ts += dTs_dt*dt


'''---------------MAIN----------------'''

'''
cur = State()

dt = 1 # 1 second
RTD = 57.2957795
dh = 15*dt/(RTD*3600)
i = 0
h0 = 1.87726272727
print "time(hr), elevation(m), Surface T(K), internal T(K), internal-ambient T (K), velocity (m/s)\n"


timerange = np.arange(start=-2.064989, stop=3.94225172727, step=dh)

time = []
elevation = []
temp_int = []
temp_dif = []
temp_surf = []
for h in timerange:
    t = (((h+h0)*RTD/15.)*3600)
    atm = fluids.atmosphere.ATMOSPHERE_1976(cur.el)
    #print atm.T
    #print colored((t/3600, cur.el, cur.Ts, cur.Ti, cur.Ti, cur.v),"yellow")
    integrate(cur,h,dt)
    if i%200==0:
        print colored((t/3600, cur.el, cur.Ts, cur.Ti, cur.Ti-atm.T, cur.v),"cyan")
        time.append(t/3600)
        elevation.append(cur.el)
        temp_int.append(cur.Ti)
        temp_dif.append(cur.Ti-atm.T)
        temp_surf.append(cur.Ts)
    if h > h0 and cur.el == 132.6:
        break
    i+=1
    #if i ==10:
    #    break

plt.figure(1)
plt.plot(time,elevation)
plt.xlabel('Time (hr)')
plt.ylabel('Elevation (m)')
plt.ylim([0,32000])

plt.figure(2)
plt.plot(time,temp_int)
plt.plot(time,temp_surf)
plt.xlabel('Time (hr)')
plt.ylabel('Temperature (K)')

plt.show()
'''


def vectorfield(w,t):
    #zdot = q1
    #T_Sdot = q2
    #T_idot = q3

    doy = 306 #temporary day of year
    lat = math.radians(35.106766) # rad

    #q1 = Elevation
    #q2 = velocity



    q1, q2, T_s, T_i, h = w
    h = math.radians(t*15)
    #if q2<132.2:
    #    q2 = 132.2
    #T_i = T_s

    bal = sphere_balloon.Sphere_Balloon(5.79,0.8)
    rad = radiation.Radiation(doy,lat,h,q1)
    q_rad  = rad.get_rad_total(lat, q1, h,5.79)
    #q_surf = bal.get_sum_q_surf(q_rad, Ts, el, v)
    #q_int  = bal.get_sum_q_int(state.Ts, state.Ti, state.el)
    #dT_sdt = (q_surf-q_int)/k


    atm = fluids.atmosphere.ATMOSPHERE_1976(q1)
    #q_int  = bal.get_sum_q_int(Ts, Ti, el)
    tm_air = atm.rho*vol*Cp_air0
    #return q_int/tm_air

    '''
    if q1 < 132.6:
        q2= 0
        '''




    f = [1,#q2,
        0,#get_acceleration(q2,q1,T_s,T_i),
        (bal.get_sum_q_surf(q_rad, T_s, q1, q2)-bal.get_sum_q_int(T_s, T_i, q1))/k,
        (bal.get_sum_q_surf(q_rad, T_s, q1, q2)-bal.get_sum_q_int(T_s, T_i, q1))/k,#bal.get_sum_q_int(T_s, T_i, q1)/tm_air,
        .45]
    return f


q1 = 132.6
q2 = 0
T_s = 287.288117979
T_i = 287.288117979
h = -2.064989

RTD = 57.2957795

abserr = 1.0e-8
relerr = 1.0e-6
stoptime= 10.0
numpoints = 250

#t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
t = np.linspace(-6, 14, 600)

#print t

w0 = [q1, q2, T_s, T_i, h]

wsol = odeint(vectorfield, w0, t,
              atol=abserr, rtol=relerr)

for t1, w1 in zip(t, wsol):
        print t1, w1[0], w1[1], w1[2], w1[3], w1[4]

print np.shape(t)
print np.shape(wsol)
#time = t[ : , 0]
T_SURFACE = wsol[ : , 2]
T_INTERNAL = wsol[ : , 3]
plt.figure(2)
plt.plot(t,T_SURFACE)
plt.plot(t,T_INTERNAL)
plt.xlabel('Time (hr)')
plt.ylabel('T)s (m)')

print "k", k
plt.show()

'''
s = State()
der = StateDerrivative()
h0 = -2.064989 #before sunrise
get_dTs(s,h0)
get_dTi(s)

eval(s,1,0,der)
print .4714*math.sqrt(0.0)*pow(0.710553430019, (1./3.))
integrate(s,-1,1)

print dh
print timerange.size
'''


'''
def firstorder(y,t,k,u):
    tau = 5.0
    dydt = (-y+k*u)/tau
    return dydt

s = State()
get_acceleration(s)
#get_dTs(s,0)
get_dTi(s)

h= np.linspace(0,math.pi,100)
#T_s = odeint(get_dTs,0,h,args=(s,))


t = np.linspace(0,10,11)
k = 2.0
u = 1.0
y = odeint(firstorder,0,t,args=(k,u))
print y
'''
