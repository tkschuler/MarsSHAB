import math
import fluids
import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pow, radians, degrees, fabs, exp, acos, asin

import config

#MARS RADIATION

class MarsRadiation:
    e = 0.09341233 #Mars Eccentiricty
    a = 1.524

    Ls = radians(config.mars_properties['Ls'])
    optical_depth = config.mars_properties['optical_depth']
    emissGround = config.mars_properties['emissGround']
    albedo = config.mars_properties['albedo']

    P0 = 669.0 #Pressure @ Surface Level (Pa)
    SB = 5.670373E-8 #Stefan_Boltzan Constant
    RE = 3376000. # Radius of Mars (m)

    emissEnv = config.balloon_properties['emissEnv']
    absEnv = config.balloon_properties['absEnv']
    #transEnv = .1
    #refEnv = .1 #revlectivity of envelope

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
        w = radians(15*h-180.)
        decl = self.get_declination(Ls)
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
        if G_bh < 0:
            G_bh = 0
        return G_bh

    def get_Mars_IR(self,el):
        '''FIX THIS'''
        T_surface = 225
        p = self.get_P(el)
        IR_trans = 1.716-0.5*(math.exp(-0.65*p/MarsRadiation.P0) + math.exp(-0.095*p/MarsRadiation.P0))
        IR_tot = IR_trans*MarsRadiation.emissGround*MarsRadiation.SB*pow(T_surface,4)
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
        power_mars_IR = mars_IR*vf*surfArea*MarsRadiation.emissEnv  #emissIR = absIR

        '''
        sky_IR = self.get_sky_IR(el)
        power_sky_IR = sky_IR*totAbs*(1.-vf)*surfArea
        '''

        rad_tot = power_direct+power_mars_IR + power_mars_IR
        return rad_tot
