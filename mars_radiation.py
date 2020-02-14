import math
import fluids
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pow, radians, degrees, fabs, exp, acos, asin

#MARS RADIATION

class MarsRadiation:
    e = 0.09341233 #Mars Eccentiricty
    Ls = 0 # Day, Vernal Equinox
    a = 1.524
    optical_depth = .5 #typical on clear days #assumption
    P0 = 669.0 #Pressure @ Surface Level (Pa)
    SB = 5.670373E-8 #Stefan_Boltzan Constant
    RE = 3376000. # Radius of Mars (m)

    emissGround = .95 #assumption
    albedo = 0.17 #assumption

    absEnv = .6 #absorbiviy of envelope
    emmisEnv = .03 #emissivity of envelope
    transEnv = .1
    refEnv = .1 #revlectivity of envelope

    '''Presure and Temperature model as function of elevation from
        https://www.grc.nasa.gov/www/k-12/airplane/atmosmrm.html'''

    def get_P(self,el):
        """Pressure at Elevation

        :param el: Elevation (m)
        :type el: float
        :returns: Pressure (Pa)
        :rtype: float
        """
        try:
            p = .699*exp(-.00009*el)*1000.
        except:
            p = 699.
        return p #pascals

    def get_g(self,el):
        """Gravity at Elevation

        :param el: Elevation (m)
        :type el: float
        :returns: Accelertation (m/s^2)
        :rtype: float
        """
        #Proportional to radius
        g= 3.711/pow(((MarsRadiation.RE+el)/MarsRadiation.RE),2)
        return g

    def get_T(self, el):
        """Temperature at Elevation

        :param el: Elevation (m)
        :type el: float
        :returns: Temperature (K)
        :rtype: float
        """
        if el < 7000:
            T = -31-.000998*el
        else:
            T = -23.4-.00222*el
        T += T + 273.15 #convert to Kelvin
        return T

    def get_rho(self,el):
        """Density at Elevation

        :param el: Elevation (m)
        :type el: float
        :returns: Density (kg/m^3)
        :rtype: float
        """
        T = self.get_T(el) - 273.5
        P = self.get_P(el) / 2000
        rho = P/(.1921*(T+273.1)) #kg/m^3
        return rho

    def get_Beam_I0(self,Ls):
        '''
        f = MarsRadiation.Ls-248 #degrees
        r = MarsRadiation.a*(1-math.pow(MarsRadiation.e,2))/(1+MarsRadiation.e*math.cos(f))
        S = 1371. # W/m^2 solar Constant
        return S/math.pow(r,2)
        '''
        #G_ob
        return 590*math.pow(1+MarsRadiation.e*cos(Ls-radians(248)),2)/pow((1-pow(MarsRadiation.e,2)),2)

    def get_declination(self,Ls):
        """Expression from http://large.stanford.edu/courses/2017/ph240/black1/docs/nasa-tm-102299.pdf

        :returns: Approximate solar declination (rad)
        :rtype: float
        """

        decl = asin(sin(radians(24.936))*sin(Ls))
        return decl

    def get_zenith(self,Ls,lat,h):
        w = math.pi/180*(15*h-180.)#radians(15*h-180.)
        #print 'w', w
        decl = self.get_declination(Ls)
        #print sin(lat)
        zen = acos(sin(lat)*sin(decl)+cos(lat)*cos(decl)*cos(w))
        return zen

    def get_air_mass(zen):
        return 1/cos(zen)

        '''
    def get_optical_depth(self,Ls,lat):
        print "asldkals", -1*pow(Ls-215,2)/730
        print 0.779*exp(pow(-(Ls-215),2)/730.)#+exp(pow(-(Ls-295.),2)/730.)
        return max([.5,16787.*(1+lat/150.)/(1917+pow((lat+38.27),2))])
        '''

    def get_surface_radiation(self,Ls, zen):
        air_mass = 1/cos(zen)
        return self.get_Beam_I0(Ls)*exp(-MarsRadiation.optical_depth*air_mass)

    def get_global_irradiance(self,Ls,lat,h):
        G_ob = self.get_Beam_I0(Ls)
        G_ob = 590
        zen = self.get_zenith(Ls,lat,h)
        G_h = G_ob*cos(zen)*.8/0.9#MarsRadiation.optical_depth/0.9
        if G_h < 0 :
            G_h = 0
        return G_h

    def get_beam_irradiance(self,Ls,lat,h):
        G_ob = 590
        zen = self.get_zenith(Ls,lat,h)
        G_bh = G_ob*cos(zen)*exp(-1*MarsRadiation.optical_depth/cos(zen))
        #print G_bh
        if G_bh < 0:
            G_bh = 0
        return G_bh

    def get_Mars_IR(self,el):
        '''FIX THIS'''
        T_surface = 225
        p = self.get_P(el)
        #print p
        IR_trans = 1.716-0.5*(math.exp(-0.65*p/MarsRadiation.P0) + math.exp(-0.095*p/MarsRadiation.P0))
        #print IR_trans
        IR_tot = IR_trans*MarsRadiation.emissGround*MarsRadiation.SB*pow(T_surface,4)
        #print "IR,tot", IR_tot
        return IR_tot

    def get_albedo_flux(self,Ls,lat,h):
        zen = self.get_zenith(Ls,lat,h)
        I_sun = self.get_global_irradiance(Ls,lat,h)
        albedo_flux = MarsRadiation.albedo*I_sun*sin(zen)
        #print "albedo flux", albedo_flux
        return albedo_flux



    def get_rad_total(self,lat,Ls,el,t,d):
        """Total Radiation as a function of elevation, time of day, and balloon surface area

        :param el: Elevation (m)
        :type el: float
        :returns: Total radiation (W/m^2)
        :rtype: float
        """

        h = t/3600.

        projArea = 0.25*math.pi*d*d
        surfArea = math.pi*d*d

        try:
            hca = math.asin(MarsRadiation.RE/(MarsRadiation.RE+el)) #half cone angle
        except:
            hca = math.asin(radians(1))

        vf = 0.5*(1. - math.cos(hca)) #viewfactor

        G_h = self.get_global_irradiance(Ls,lat,h)
        power_direct = G_h*MarsRadiation.absEnv*projArea

        '''
        diffuse_I = self.get_diffuse_SI(zen, el)
        power_diffuse = diffuse_I*totAbs*(1.-vf)*surfArea
        '''

        albedo_flux = self.get_albedo_flux(Ls,lat,h)
        power_reflected = albedo_flux*MarsRadiation.absEnv*vf*surfArea

        mars_IR = self.get_Mars_IR(el)
        power_mars_IR = mars_IR*vf*surfArea*MarsRadiation.emmisEnv  #emissIR = absIR

        '''
        sky_IR = self.get_sky_IR(el)
        power_sky_IR = sky_IR*totAbs*(1.-vf)*surfArea
        '''

        rad_tot = power_direct+power_mars_IR + power_mars_IR
        return rad_tot

m = MarsRadiation()

Ls = radians(153)
lat = radians(22.3)
h = 12.
d = 15 #m
'''
print m.get_global_irradiance(Ls,lat,h)

h_ang = np.arange(start=5, stop=19, step=.1)

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
plt.show()
'''

'''
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
        #zen = r.get_zenith(lat,math.radians(h*15))
        #I = r.get_direct_SI(zen,el) + r.get_diffuse_SI(zen,el) + r.get_reflected_SI(zen,el)
        #IR = r.get_earth_IR(el)+r.get_sky_IR(el)
        #total_radiation = I+IR


        total_radiation = m.get_global_irradiance(Ls,lat,h) + m.get_albedo_flux(Ls,lat,h) + m.get_Mars_IR(el)
        #total_radiation = m.get_rad_total(lat,Ls,el,h,d)/100

        #el = y[j]/1000 #convert to km

        print h, el, total_radiation
        z[j,i] = total_radiation

for i in range(0,x.size):
    h = x[i]
    direct.append(m.get_global_irradiance(Ls,lat,h))
    albedo.append(m.get_albedo_flux(Ls,lat,h))
    IR.append(m.get_Mars_IR(0))

z_min, z_max = np.amin(z), np.amax(z)
#z_min, z_max = 0, 600

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

plt.figure(2)
plt.title('Total Solar Intensity at Martian Surface Level')
plt.plot(x,direct, label="Direct Solar Flux")
plt.plot(x,albedo, label="Reflected Albedo Flux")
plt.plot(x,IR, label="IR Radiation from Surface")
plt.xlabel('Time (hr)')
plt.ylabel('Intensity (W/m^2)')
plt.legend(loc='upper right')
plt.show()
'''

'''
el = np.arange(start=0, stop=30000, step=5)
T = []
for e in np.nditer(el):
    T.append(m.get_T(e))
plt.figure(1)
plt.plot(T,el)
plt.xlabel('Temp (K)')
plt.ylabel('Elevation (m)')
plt.show()
'''
