import math
import fluids
import numpy as np
import matplotlib.pyplot as plt

#EARTH RADIATION

class Radiation:
    #Radiation Constants
    I0 = 1358 #Direct Solar Radiation Level
    e = 0.016708 #Eccentricity of Earth's Orbit
    P0 = 101325 #Standard Atmospheric Pressure at Sea Level
    cloudElev = 3000 #m
    cloudFrac = 0.0 #percent cloud coverage [0,1]
    cloudAlbedo = .6
    albedoGround = .2 #ground albedo [0,1]
    tGround = 293
    emissGround = .95
    SB = 5.670373E-8
    RE = 6371000 #m Radius of Earth
    radAbs = .8
    emissEnv = radAbs
    radRef= .1
    radTrans = .1

    def __init__(self, doy, lat, h_ang, el):
        self.doy = doy
        self.lat = lat
        self.h_ang = h_ang
        self.el = el

    def get_SI0(self):
        """ Incident solar radiation

        :returns: The incident solar radiation above Earths atm (W/m^2)
        :rtype: float
        """
        f = 2*math.pi*self.doy/365 #true anomaly
        e2 = pow(((1.+Radiation.e)/(1.-Radiation.e)),2) -1.
        return Radiation.I0*(1.+0.5*e2*math.cos(f))

    def get_declination(self):
        """Expression from http://en.wikipedia.org/wiki/Position_of_the_Sun

        :returns: Approximate solar declination (rad)
        :rtype: float
        """
        return -.4091*math.cos(2*math.pi*(self.doy+10)/365)

    def get_zenith(self,lat, h_ang):
        """get zenith angle

        :param lat: Lattitude (rad)
        :type lat: float
        :param h_ang: Solar Hour Angle (rad)
        :type h_ang: float
        :returns: The approximate solar hour angle
        :rtype: float
        """
        decl = self.get_declination()
        return math.acos(math.sin(self.lat)*math.sin(decl)+math.cos(self.lat)*math.cos(decl)*math.cos(h_ang))


    def get_air_mass(self,zen, el):
        """Air Mass at elevation

        :param zen: Solar Angle (rad)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :returns: The approximate air mass (unitless)
        :rtype: float
        """

        atm = fluids.atmosphere.ATMOSPHERE_1976(el)
        p = atm.P #pressure at current elevation
        am = (p/Radiation.P0)*(math.sqrt(1229 + pow((614*math.cos(zen)),2))-614*math.cos(zen))
        return am

    def get_trans_atm(self,zen,el):
        """get zenith angle

        :param zen: Solar Angle (rad)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :returns: The atmospheric trasmittance (unitless)
        :rtype: float
        """

        if math.fabs(zen) > math.pi/2.:
            return 0.0
    	am = self.get_air_mass(zen, el)
        return 0.5*(math.exp(-0.65*am) + math.exp(-0.095*am))

    def get_direct_SI(self,zen,el):
        """Get Direct Solar Radiation

        :param zen: Solar Angle (rad)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :returns: Tntensity of the direct solar radiation (W/m^2)
        :rtype: float
        """

        SI0 = self.get_SI0()
        trans = self.get_trans_atm(zen, el)
        return trans*SI0

    def get_diffuse_SI(self,zen,el):
        """Diffuse Solar Radiation from sky

        :param zen: Solar Angle (rad)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :returns: The intensity of the diffuse solar radiation from the sky (W/m^2)
        :rtype: float
        """

        if(zen > math.pi/2.):
            return 0.0
        SI0 = self.get_SI0()
        trans = self.get_trans_atm(zen, el)
        if el < Radiation.cloudElev:
            return (1-Radiation.cloudFrac)*0.5*SI0*math.sin(math.pi/2.-zen)*(1.-trans)/(1-1.4*math.log(trans))
        else:
            return 0.5*SI0*math.sin(math.pi/2.-zen)*(1.-trans)/(1-1.4*math.log(trans))


    def get_reflected_SI(self,zen,el):
        """Diffuse Solar Radiation from sky

        :param zen: Solar Angle (rad)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :returns: The intensity solar radiation reflected by the Earth (W/m^2)
        :rtype: float
        """

        if(zen > math.pi/2.):
             return 0.0
        incident_SI = self.get_SI0()
        tau_atm = self.get_trans_atm(zen,el)
        if el < Radiation.cloudElev:
            albedo = (1.-Radiation.cloudFrac)*Radiation.albedoGround;
        else:
            albedo = (1.-Radiation.cloudFrac)*(1-Radiation.cloudFrac)*Radiation.albedoGround + Radiation.cloudAlbedo*Radiation.cloudFrac
        return albedo*tau_atm*incident_SI*math.sin(math.pi/2.-zen)

    def get_earth_IR(self,el):
        """Infared Radiation from Earth's surface

        :param el: Elevation (m)
        :type el: float
        :returns: Intensity of IR radiation emitted from earth (W/m^2)
        :rtype: float
        """
        atm = fluids.atmosphere.ATMOSPHERE_1976(el)
        p = atm.P #pressure at current elevation
        IR_trans = 1.716-0.5*(math.exp(-0.65*p/Radiation.P0) + math.exp(-0.095*p/Radiation.P0))
        if el < Radiation.cloudElev:
            tEarth = Radiation.tGround
        else:
            clouds = fluids.atmosphere.ATMOSPHERE_1976(Radiation.cloudElev)
            tEarth = Radiation.tGround*(1.-Radiation.cloudFrac) + clouds.T*Radiation.cloudFrac
        return IR_trans*Radiation.emissGround*Radiation.SB*pow(tEarth,4)

    def get_sky_IR(self,el):
        """Infared Radiation from Sky

        :param el: Elevation (m)
        :type el: float
        :returns: Intensity of IR radiation emitted from sky (W/m^2)
        :rtype: float
        """
        return np.fmax(-0.03*el+300.,50.0)

    def get_rad_total(self,lat,el,h,d):
        """Total Radiation as a function of elevation, time of day, and balloon surface area

        :param el: Elevation (m)
        :type el: float
        :returns: Total radiation (W/m^2)
        :rtype: float
        """

        #some constant things
        #this doesn't make sense
        radRef = Radiation.radRef + Radiation.radRef*Radiation.radRef +  Radiation.radRef*Radiation.radRef*Radiation.radRef
        totAbs = Radiation.radAbs + Radiation.radAbs*Radiation.radTrans + Radiation.radAbs*Radiation.radTrans*radRef
        projArea = 0.25*math.pi*d*d
        surfArea = math.pi*d*d
        #--------------------------------------------------------------------------

        hca = math.asin(Radiation.RE/(Radiation.RE+el)) #half cone angle
        #print "el", el
        #print "hca: ", hca
        vf = 0.5*(1. - math.cos(hca)) #viewfactor

        zen = self.get_zenith(self.lat, h)
        direct_I = self.get_direct_SI(zen, el)
        power_direct = direct_I*totAbs*projArea

        diffuse_I = self.get_diffuse_SI(zen, el)
        power_diffuse = diffuse_I*totAbs*(1.-vf)*surfArea

        reflected_I = self.get_reflected_SI(zen, el)
        power_reflected = reflected_I*totAbs*vf*surfArea

        earth_IR = self.get_earth_IR(el)
        power_earth_IR = earth_IR*totAbs*vf*surfArea

        #print "totAbs", totAbs, "vf," , vf, "surfArea", surfArea

        sky_IR = self.get_sky_IR(el)
        power_sky_IR = sky_IR*totAbs*(1.-vf)*surfArea

        #print "\npower_sky_IR", power_sky_IR, "power_earth_IR", power_earth_IR, "power_reflected", power_reflected, "power diffuse", power_diffuse, "power_direct", power_direct

        rad_tot = (power_direct+power_diffuse+power_reflected+power_earth_IR+power_sky_IR)  #Somewhere along the way this became a factor of 100 greater?
        #print rad_tot
        return rad_tot


doy = 306 #temporary day of year
lat = math.radians(35.106766) # rad
h_ang = 0
el = 0 #elevation (m)

r = Radiation(doy,lat,-2.064989,136.6)
q_rad  = r.get_rad_total(lat, 136.6, -2.064989,5.79)
