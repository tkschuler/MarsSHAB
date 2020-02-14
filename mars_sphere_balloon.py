import math
import fluids
import numpy as np
import matplotlib.pyplot as plt
import mars_radiation
from termcolor import colored
from math import radians

from scipy.integrate import odeint
from scipy import optimize
from math import pow, fabs

class Mars_Sphere_Balloon:
    Cp_co2 = 735.0 #J/Kg*K #assumpiton?
    Rsp_co2 = 44.010
    RE = 3376000.0 #m Radius of Mars
    SB = 5.670373E-8
    #M =  11.3 #kg
    #cf = 910 #specific heat of aluminim
    cf = 320.0 #specific heat of germanium
    cv_Co2 = 735.0 #J/Kg*K #assumpiton?

    def __init__(self, d, emissEnv):
        self.d = d #diameter
        self.surfArea = math.pi*d*d
        self.emissEnv = emissEnv
        self.massEnv = self.surfArea*.009 #density of material

        self.vol = math.pi*4/3*pow((d/2),3)

    def get_dynamic_viscocity_co2(self,T):
        """Returns Dynamic Viscocity of CO2 as function of Temperature
        source: https://rotorcraft.arc.nasa.gov/Publications/files/koning_AIAA_2019.pdf
        source: https://www.lmnoeng.com/Flow/GasViscosity.php

        :param el: Temperaure (K)
        :type el: float
        :returns: Dynamic Viscocity ()
        :rtype: float
        """

        C = 240
        To = 527.67
        T = T*(9/5)
        a = 0.555*To+C
        b = 0.555*T + C
        mu = 0.01480*(a/b)*math.pow((T/To),3/2)/1000
        '''CONSTANT DUE TO LOW PRESSURE/TEMPERATURE'''
        return 1.130E-5 #interpolated

    def get_k_co2(self,T):
        """Returns Thermal Conductivity of CO2 as function of Temperature
        source: https://journals.ametsoc.org/doi/pdf/10.1175/BAMS-D-12-00158.1

        :param el: Temperaure (K)
        :type el: float
        :returns: Thermal Conductivity ()
        :rtype: float
        """

        #k = 1.1691*11.9E-40E-7*pow(T,2)+1.3327*10E-5*T+2.2469*10E-3
        #print "T", T
        #k = 0.0241*math.pow((T/273.15),0.9)
        #https://www.engineeringtoolbox.com/carbon-dioxide-thermal-conductivity-temperature-pressure-d_2019.html
        '''CONSTANT DUE TO LOW PRESSURE/TEMPERATURE'''
        k = 10.024/1000 #interpolated
        return k


    def get_Pr_co2(self,T):
        """Returns Prandtl Number of CO2 as function of Temperature
        source: https://www.engineeringtoolbox.com/carbon-dioxide-prandtl-number-viscosity-heat-capacity-thermal-conductivity-d_2024.html

        :param el: Temperaure (K)
        :type el: float
        :returns: Thermal Conductivity ()
        :rtype: float
        """

        '''CONSTANT DUE TO LOW PRESSURE/TEMPERATURE'''
        Pr = 0.76 #interpolated
        return Pr

        '''------------------------SOLVE FOR SURFACE TEMPERATURE----------------------------------------'''


    def get_Nu_ext(self,Ra, Re, Pr):
        """External Nusselt Number

        :param Ra: Raleigh's number
        :type Ra: float
        :param Re: Reynold's number
        :type Re: float
        :param Pr: Prandtl Number
        :type Pr: float
        :returns: External Nusselt Number
        :rtype: float
        """

        Nu_n = 0.0
        if Ra < 1.5E8:
            try:
                Nu_n = 2.0 + 0.6*pow(Ra,0.25)
            except:
                Nu_f = 2
        else:
            Nu_n = 0.1*pow(Ra, 0.34)
        Nu_f = 0.0
        if Re < 5E4:
            #print colored(("Ra: ", Ra, "Re:", Re, "Pr: ", Pr), "red")
            try:
                Nu_f = 2 + 0.47*math.sqrt(Re)*pow(Pr, (1./3.))
            except:
                Nu_f = 2
        else:
            Nu_f = (0.0262*pow(Re, 0.8) - 615.)*pow(Pr, (1./3.));
        return np.fmax(Nu_f, Nu_n);

    def get_Nu_free(self,T,el,Pr):
        m = mars_radiation.MarsRadiation()
        T_atm = m.get_T(el)
        p_atm = m.get_P(el)
    	rho_atm = m.get_rho(el)
        g = m.get_g(el)

        Pr_atm = self.get_Pr_co2(T_atm)

    	T_avg = 0.5*(T_atm + T)
        rho_avg = p_atm/(Mars_Sphere_Balloon.Rsp_co2*T_avg)
        Pr_avg = self.get_Pr_co2(T_avg)

        exp_coeff = 1./T_avg
        mu = self.get_dynamic_viscocity_co2(T_avg)

        Gr = (pow(rho_atm,2)*g*fabs(T-T_atm)*pow(self.d,3))/(T_atm*pow(mu,2))

        #print "Gr,", Gr, "Pr," , Pr
        try:
            Nu = 2 + .45*pow((Gr*Pr),.25)
        except:
            Nu = 2.45
        return Nu

        '''
        # This was from the Bovine Paper, not using it.
    def get_Nu_forced(self,Re, Ra, Pr):
        Nu_n = 0.0
        if Ra < 1.5E8:
            Nu_n = 2.0 + 0.6*pow(Ra,0.25)
        else:
            Nu_n = 0.1*pow(Ra, 0.34)
        Nu_f = 0.0
        if Re < 5E4:
            #print colored(("Ra: ", Ra, "Re:", Re, "Pr: ", Pr), "red")
            try:
                Nu_f = 2 + 0.47*math.sqrt(Re)*pow(Pr, (1./3.))
                #print colored(Nu_f, "yellow")
            except:
                Nu_f = 2
                #print colored("WTFFFFFFF why is there a math domain error", "yellow")
        else:
            Nu_f = (0.0262*pow(Re, 0.8) - 615.)*pow(Pr, (1./3.));
        return np.fmax(Nu_f, Nu_n);
        '''

    def get_q_ext(self, T_s, el, v):
        """External Heat Transfer

        :param zen: Surface Temperature of Envelope (K)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :param el: velocity (m/s)
        :type el: float
        :returns: Power transferred from sphere to surrounding atmosphere due to convection(W)
        :rtype: float
        """

        m = mars_radiation.MarsRadiation()
        T_atm = m.get_T(el)
        p_atm = m.get_P(el)
    	rho_atm = m.get_rho(el)
        g = m.get_g(el)

    	Pr_atm = self.get_Pr_co2(T_atm)

    	T_avg = 0.5*(T_atm + T_s)
        rho_avg = p_atm/(Mars_Sphere_Balloon.Rsp_co2*T_avg)
        Pr_avg = self.get_Pr_co2(T_avg)

        exp_coeff = 1./T_avg
        kin_visc = self.get_dynamic_viscocity_co2(T_avg)/rho_avg

        k = self.get_k_co2(T_avg)
        alpha = k/(rho_avg*Mars_Sphere_Balloon.Cp_co2)
        Ra = g*exp_coeff*pow(self.d,3)/(kin_visc*alpha)*math.fabs(T_s-T_atm)
        # Reynolds number has to be positive, therefore convert negative velocities.
    	Re = rho_atm*fabs(v)*self.d/self.get_dynamic_viscocity_co2(T_atm)

        Nu = self.get_Nu_free(T_s,el,Pr_avg)
    	k = self.get_k_co2(T_avg)

        '''External Free Convection'''
    	h = (Nu*k)/self.d
        '''External Forced Convection'''
        h_forced =  k/self.d*(2+.41*np.power(Re,0.55))
        #Take maximum value between free and forced
        h = np.fmax(h,h_forced)

        q_conv = h*self.surfArea*(T_s-T_atm)
    	return q_conv

    '''------------------------SOLVE FOR T INT----------------------------------------------'''


    def get_q_int(self,T_s, T_i, el):

        m = mars_radiation.MarsRadiation()
        T_atm = m.get_T(el)
        p_atm = m.get_P(el)
    	rho_atm = m.get_rho(el)
        g = m.get_g(el)

    	Pr = self.get_Pr_co2(T_i)


    	mu = self.get_dynamic_viscocity_co2(T_i)
    	k = self.get_k_co2(T_i)
    	h = 0.13*k*pow((pow(rho_atm,2)*g*fabs(T_s-T_i)*Pr)/(T_i*pow(mu,2)),(1/3))
        #h = .01
        #print "h", h
        q_int = h*self.surfArea*(T_s-T_i)
        #print colored(("q_int", q_int),"yellow")
    	return q_int

    def get_sum_q_surf(self,q_rad, T_s,T_i, el, v):
        """External Heat Transfer

        :param q_rad: Power input from external radiation (W)
        :type q_rad: float
        :param T_s: Surface Temperature of Envelope (K)
        :type T_s: float
        :param el: Elevation (m)
        :type el: float
        :param v: velocity (m/s)
        :type v: float
        :returns: The sum of power input to the balloon surface (W)
        :rtype: float
        """
        # https://www.sciencedirect.com/science/article/pii/S0038092X15002418
        # http://www.ae.utexas.edu/courses/ase261/balloon/BalloonTrajectory.pdf

        q_ce = -self.get_q_ext(T_s, el, v) #Heat Loss due to External Convection
    	q_re = -self.emissEnv*Mars_Sphere_Balloon.SB*self.surfArea*(pow(T_s,4)) #Heat Loss due to radiation
        q_ci = -self.get_q_int(T_s, T_i, el) #Heat Transfer due to Internal Convection
    	return q_rad + q_ce + q_re + q_ci

    def solve_T_surf(self,q_rad, el, v):
        def f(T_s):
            return q_rad-self.emissEnv*Mars_Sphere_Balloon.SB*pow(T_s,4)*self.surfArea
        T_s = optimize.newton(f, 400)
        if T_s< 211.5:
            T_s = 211.5
        return T_s

    '''

    def get_Nu_int(self,Ra):
        print "RA", Ra
        try:
        	if Ra < 1.35E8:
        		return 2.5*(2+0.6*pow(Ra,0.25))
        	else:
        		return 0.325*pow(Ra, 0.333)
        except:
            print colored("negative exponent", "red")
            return 0.0

    def get_q_int(self,T_s, T_i, el):

        m = mars_radiation.MarsRadiation()
        T_atm = m.get_T(el)
        p_atm = m.get_P(el)
    	rho_atm = m.get_rho(el)
        g = m.get_g(el)

    	T_avg = 0.5*(T_s+T_i)
    	rho_avg = p_atm/(Mars_Sphere_Balloon.Rsp_co2*T_avg)
    	Pr = self.get_Pr_co2(T_avg)
    	exp_coeff = 1./T_avg
    	kin_visc = self.get_dynamic_viscocity_co2(T_avg)/rho_avg
    	Ra = self.get_Pr_co2(T_atm)*g*math.fabs(T_i-T_s)*pow(self.d,3)*exp_coeff/(kin_visc*kin_visc)
    	Nu = self.get_Nu_int(Ra)
    	k = self.get_k_co2(T_avg)
    	h = Nu*k/self.d
    	return h*self.surfArea*(T_i-T_s)



    def get_sum_q_int(self, T_s, T_i, el):
        q_ci = -self.get_q_int(T_s, T_i, el)
        #should there even be IR transfer between internal air & balloon surf?
        #i dont think so
        #double q_ri = 0*E_int*SB_CONST*(pow(T_s,4)-pow(T_i,4))*surface_area;
        return q_ci #// + q_ri;

'''

Ls = radians(153)
lat = radians(22.3)
#h = 5
d = 18 #m
el = 0.
v = 0

m = Mars_Sphere_Balloon(d,.03)
r = mars_radiation.MarsRadiation()

'''
T_atm = 211.5
r = mars_radiation.MarsRadiation()
h_ang = np.arange(start=0, stop=20, step=.1)
rad_tot = []
T_s = []
for h in np.nditer(h_ang):
    #h = math.radians(h)
    rad_total = r.get_rad_total(lat,Ls,el,h,d)
    print "RAD_TOTAL", rad_total
    rad_tot.append(rad_total)
    T_s.append(m.solve_T_surf(rad_total,el,v))
    #T_i.append(m.solve_T_i(T_s, el))

plt.figure(1)
plt.plot(h_ang,T_s)
print "CONVECTION:", m.get_external_free_convection(200,0)
plt.show()
'''

'''
def model(T_s,t):
    #h = t/3600.
    r = mars_radiation.MarsRadiation()
    rad_total = r.get_rad_total(lat,Ls,el,t,d)
    dTdt = m.get_sum_q_surf(rad_total, T_s, el, v)/((m.cf*m.massEnv))
    return dTdt

# initial condition
T_s0 = 200
# time points
t = np.linspace(0,100*3600,100*3600)
# solve ODE
T_s = odeint(model,T_s0,t)
'''
'''
dt = 1
T_s = [211.5]
T_i = [212.3]
el = 100.0
v = 0.0

t = np.linspace(0,100*3600,100*3600)
for i in range(0,len(t)-1):
    T_s_cur= T_s[i]
    T_i_cur = T_i[i]

    print T_s_cur, T_i_cur
    r = mars_radiation.MarsRadiation()
    rad_total = r.get_rad_total(lat,Ls,el,t[i],d)
    dT_sdt = m.get_sum_q_surf(rad_total, T_s_cur, T_i_cur, el, v)/((m.cf*m.massEnv))
    rho_int = r.get_rho(T_i_cur)#r.get_rho(T_i[i])

    tm_air = rho_int*m.vol*m.Cp_co2
    dT_idt = m.get_q_int(T_s_cur, T_i_cur, el)/(tm_air)

    T_s.append(T_s_cur+dt*dT_sdt)
    T_i.append(T_i_cur+dt*dT_idt)




# plot results
plt.plot(t/3600,T_s)
plt.plot(t/3600,T_i)
plt.xlabel('Time (hr)')
plt.ylabel('T_s (K)')
#print "k", ((m.cf*m.M))

print r.get_rad_total(lat,Ls,el,100*3600,d)
plt.show()
'''
