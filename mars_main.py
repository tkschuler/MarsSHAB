import math
import fluids
import numpy as np
import mars_radiation
import mars_sphere_balloon
from termcolor import colored
from math import radians
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from sympy import symbols, diff, sin

import config

class Mars_Main:

    def __init__(self):
        self.Cp_co2 = config.mars_properties['Cp_co2']
        self.Cv_co2 = config.mars_properties['Cv_co2']
        self.Rsp_co2 = config.mars_properties['Rsp_co2']
        self.Ls = radians(config.mars_properties['Ls'])
        self.lat = radians(config.mars_properties['lat'])

        self.d = config.balloon_properties['d']
        self.emissEnv = config.balloon_properties['emissEnv']
        self.cp = config.balloon_properties['cp']

        self.vol = math.pi*4/3*pow((self.d/2),3) #volume m^3
        self.surfArea = math.pi*self.d*self.d #m^2
        self.cs_area = math.pi*self.d*self.d/4.0 #m^2

        self.vm_coeff = .5 #virtual mass coefficient
        self.massEnv = self.surfArea*(9./1000.) #kg
        self.mdot = 0 # Initial Mass Flow rate
        self.k = self.massEnv*self.cp #thermal mass coefficient

        self.dt = config.dt

    def get_acceleration(self,v,el,T_s,T_i,mp):
        m = mars_radiation.MarsRadiation()
        T_atm = m.get_T(el)
        p_atm = m.get_P(el)
        rho_atm = p_atm/(189*T_atm) #m.get_rho(el)
        g = m.get_g(el)

        rho_int = p_atm/(self.Rsp_co2*T_i)
        F_b = (rho_atm - rho_int)*self.vol*g
        Cd = .8
        F_d =  Cd*(0.5*rho_atm*math.fabs(v)*v)*self.cs_area
        vm = (self.massEnv + mp) + rho_atm*self.vol + self.vm_coeff*rho_atm*self.vol
        accel = (F_b  - F_d - (self.massEnv + mp)*g)/vm
        return accel

    def get_convection_vent(self,T_i,el):
        m = mars_radiation.MarsRadiation()
        T_atm = m.get_T(el)
        Q_vent =  self.mdot*self.Cv_co2*(T_i-T_atm)
        return Q_vent

    def solve_states(self,T_s,T_i,el,v,mp,vent,timespace,alt_dz,vel_dz):
        for i in range(0,len(timespace)-1):
            bal = mars_sphere_balloon.Mars_Sphere_Balloon()
            rad = mars_radiation.MarsRadiation()
            q_rad  = rad.get_rad_total(self.lat,self.Ls,el[i],timespace[i],self.d)

            rho_atm = rad.get_rho(el[i])

            p_atm = rad.get_P(el[i])
            rho_int = p_atm/(self.Rsp_co2*T_i[i])
            tm_air = rho_int*self.vol*self.Cp_co2

            dT_sdt = (bal.get_sum_q_surf(q_rad, T_s[i], T_i[i], el[i], v[i]))/(bal.cf*self.massEnv)
            dT_idt = (bal.get_q_int(T_s[i], T_i[i], el[i])-self.get_convection_vent(T_i[i],el[i]))/tm_air
            T_s.append(T_s[i]+dT_sdt*self.dt)
            T_i.append(T_i[i]+dT_idt*self.dt)

            dzdotdt = self.get_acceleration(v[i],el[i],T_s[i],T_i[i],mp)
            zdot = v[i] + dzdotdt*self.dt
            z = el[i]+zdot*self.dt#+dzdotdt*pow(dt,2)#*1000

            if z<0:
                v.append(0)
                el.append(0)
            else:
                v.append(zdot)
                el.append(z)

            if el[i] > alt_dz:
                self.mdot = vent

            if v[i] < vel_dz:
                self.mdot = 0

            if i%3600 == 0:
                print "t", timespace[i+1]/3600,  "el", el[i+1], "v", v[i+1], "accel", dzdotdt, "T_s", T_s[i+1], "T_i", T_i[i+1]

        return [T_s,T_i,el,v]
