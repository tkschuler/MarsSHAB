import math
import fluids
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import radiation
from termcolor import colored

class Sphere_Balloon:
    Cp_air0 = 1003.8
    Rsp_air = 287.058
    SB_CONST =  5.670373E-8

    def __init__(self, d, emissEnv):
        self.d = d #diameter
        self.surfArea = math.pi*d*d
        self.emissEnv = emissEnv

    def get_viscocity(self,T):
        return 1.458E-6*math.pow(T,1.5)/(T+110.4)

    def get_conduction(self,T):
        return 0.0241*math.pow((T/273.15),0.9)

    def get_Pr(self,T):
        k = self.get_conduction(T) #Thermal diffusivity
        return self.get_viscocity(T)*Sphere_Balloon.Cp_air0/k

    '''------------------------SOLVE FOR T_S----------------------------------------------'''

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

        atm = fluids.atmosphere.ATMOSPHERE_1976(el)

        T_atm = atm.T
        p_atm = atm.P
    	rho_atm = atm.rho
        g = atm.g

        '''	double T_atm = atm->get_T(el);
   double p_atm = atm->get_P(el);
	double rho_atm = atm->get_rho(el);
	double Pr_atm = atm->get_Pr(T_atm);

	double T_avg = 0.5*(T_atm + T_s);
   double rho_avg = p_atm/(atm->Rsp_air*T_avg);
   double Pr_avg = atm->get_Pr(T_avg);
   double exp_coeff = 1./T_avg;
   double kin_visc = atm->get_visc(T_avg)/rho_avg;
   double Ra = Pr_avg*atm->get_g(el)*fabs(T_s-T_atm)*pow(diameter,3)*exp_coeff/(kin_visc*kin_visc);'''

    	Pr_atm = self.get_Pr(T_atm)

    	T_avg = 0.5*(T_atm + T_s)
        rho_avg = p_atm/(Sphere_Balloon.Rsp_air*T_avg)
        Pr_avg = self.get_Pr(T_avg)

        exp_coeff = 1./T_avg;
        kin_visc = self.get_viscocity(T_avg)/rho_avg

        #Not sure if Raleighs number is the right equation here:
        #print "Pr_avg", Pr_avg, "kin_visc", kin_visc, "exp_coeff", exp_coeff, "T_atm", T_atm, "T_s", T_s
        Ra = Pr_avg*g*math.fabs(T_s-T_atm)*math.pow(self.d,3)*exp_coeff/(kin_visc*kin_visc)
        #print "Ra", Ra
    	Re = rho_atm*v*self.d/self.get_viscocity(T_atm)
    	Nu = self.get_Nu_ext(Ra, Re, Pr_atm)
    	k = self.get_conduction(T_avg)
        #print "Ra: ", Ra, "Re:", Re, "Pr: ", Pr_avg, "v: " , v, "el: ", el
    	h = Nu*k/self.d
        #print "h:", h
        #print "Nu: ", Nu, "k: ", k, "d:", self.d
    	return h*self.surfArea*(T_s-T_atm)

    def get_sum_q_surf(self,q_rad, T_s, el, v):
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

        q_ce = -self.get_q_ext(T_s, el, v)
    	q_re = -self.emissEnv*Sphere_Balloon.SB_CONST*pow(T_s,4)*self.surfArea
        print "T_s: ", T_s, "q_rad: ", q_rad, "q_ce: ", q_ce,  "q_re: ", q_re, "v", v
        #print q_rad + q_ce + q_re
    	return q_rad + q_ce + q_re


    def solve_T_surf(self,q_rad, el, v):
        """
        :returns: The sum of power input to the balloon surface (W)
        :rtype: float
        """
    	dT = 1. #initial guesses
        atm = fluids.atmosphere.ATMOSPHERE_1976(el)
        #print "T_atm: ", atm.T
    	T_s = atm.T+10 #T_s > T_atm
    	for i in range(0,10):  #10 iterations should be enough
            #print "T_S: ", T_s, "dT: ", dT
            q2 = self.get_sum_q_surf(q_rad, T_s+dT, el, v)
            q1 = self.get_sum_q_surf(q_rad, T_s, el, v)
            dqdT = (q2-q1)/dT
            dT = q1/dqdT
            T_s -= dT
            if math.fabs(dT) < 1E-10:
                break
        return T_s


    '''------------------------SOLVE FOR T INT----------------------------------------------'''

    def get_Nu_int(sef,Ra):
	if Ra < 1.35E8:
		return 2.5*(2+0.6*pow(Ra,0.25))
	else:
		return 0.325*pow(Ra, 0.333)

    def get_q_int(self,T_s, T_i, el):

        atm = fluids.atmosphere.ATMOSPHERE_1976(el)
        T_atm = atm.T
        p_atm = atm.P

    	T_avg = 0.5*(T_s+T_i)
    	rho_avg = p_atm/(Sphere_Balloon.Rsp_air*T_avg)
    	Pr = self.get_Pr(T_avg)
    	exp_coeff = 1./T_avg
    	kin_visc = self.get_viscocity(T_avg)/rho_avg
    	Ra = self.get_Pr(T_atm)*atm.g*math.fabs(T_i-T_s)*pow(self.d,3)*exp_coeff/(kin_visc*kin_visc)
    	Nu = self.get_Nu_int(Ra)
    	k = self.get_conduction(T_avg)
    	h = Nu*k/self.d
    	return h*self.surfArea*(T_i-T_s)

    def get_sum_q_int(self, T_s, T_i, el):
        q_ci = -self.get_q_int(T_s, T_i, el)
        #should there even be IR transfer between internal air & balloon surf?
        #i dont think so
        #double q_ri = 0*E_int*SB_CONST*(pow(T_s,4)-pow(T_i,4))*surface_area;
        return q_ci #// + q_ri;


    def solve_T_i(self,T_s, el):
        T_i = T_s + 10
        dT = 1
    	for i in range(0,10):  #10 iterations should be enough
            #print "T_S: ", T_s, "dT: ", dT
            q2 = self.get_sum_q_int(T_s, T_i+dT, el)
            q1 = self.get_sum_q_int(T_s, T_i, el)
            dqdT = (q2-q1)/dT
            dT = q1/dqdT
            T_i -= dT
            if math.fabs(dT) < 1E-10:
                break
        #print T_i
        return T_i


    def get_Cd(v,el):
        return 0.8
        '''
        Re = atm->get_rho(el)*fabs(v)*diameter/atm->get_visc(atm->get_T(el));
        unsigned int i = 0;
       for(i = 0; i < cd_arr.size(); i++) {
          if(cd_arr[i].Re>Re) break;
       }
       if(i==0) return cd_arr[0].Cd;
       double Cd = (cd_arr[i].Cd - cd_arr[i-1].Cd)/(cd_arr[i].Re - cd_arr[i-1].Re);
       Cd *= (Re - cd_arr[i-1].Re);
       Cd += cd_arr[i-1].Cd;
       return Cd;
       '''


doy = 306 #temporary day of year
lat = math.radians(35.106766) # rad
h_ang = 0
el = 0 #elevation (m)

r = radiation.Radiation(doy,lat,h_ang,el)

T_s = 0
el = 0
v =0
atm = fluids.atmosphere.ATMOSPHERE_1976(el)
T_atm = atm.T

d= 5.79
emissEnv = 0.8

s = Sphere_Balloon(5.79,emissEnv)

h_ang = np.arange(start=-90, stop=90, step=.1)
h_ang2 = []
T_s = []
rad_tot = []
for h in np.nditer(h_ang):
    h = math.radians(h)
    rad_total = r.get_rad_total(lat,el,h,d)
    rad_tot.append(r.get_rad_total(lat,el,h,d))
    h_ang2.append(1./15.*math.degrees(h))
    T_s.append(s.solve_T_surf(rad_total,132,0))
    print s.solve_T_surf(rad_total,132,0)

plt.figure(1)
plt.plot(h_ang2,T_s)
plt.show()


'''
s = Sphere_Balloon(5.79,.03)
print s.get_q_ext(297,0,0)
print s.get_sum_q_surf(600, 297, 0, 0)


s = Sphere_Balloon(5.79,.03)
'''
s = Sphere_Balloon(5.79,.03)
'''
def model(T_s,t):
    r = radiation.Radiation(doy,lat,h_ang,el)
    rad_total = r.get_rad_total(lat,el,h,d)
    dTdt = s.get_sum_q_surf(rad_total, T_s, el, v)#((m.cf*m.M))
    return dTdt

T_s0 = 210
# time points
t = np.linspace(0,10)
# solve ODE
T_s = odeint(model,T_s0,t)

# plot results
plt.plot(t,T_s)
plt.xlabel('time')
plt.ylabel('y(t)')
#print "k", ((m.cf*m.M))
plt.show()
'''



#rad_total = r.get_rad_total(lat,132,-1.5,d)
#print rad_total
