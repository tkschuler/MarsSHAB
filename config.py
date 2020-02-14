
balloon_properties = dict(
    shape = 'sphere',
    d = 18, # (m) Diameter of Sphere Balloon
    mdot = 0, # (kg/s) Vent Mass Flow Rate
    cp = 320.0, #(J/(kg K)) Specific heat of envelope material

    emissEnv = .03, # Emisivity of Enevelope
    absEnv = .6, # Absorbiviy of envelope
    transEnv = .1,
    refEnv = .1 #revlectivity of envelope
)

mars_properties = dict(
    Cp_co2 = 735.0, # (J/Kg*K)  Specifc Heat Capacity, Constant Pressure
    Cv_co2 = 657.0, # (J/Kg*K)  Specifc Heat Capacity, Constant Volume
    Rsp_co2 = 188.92, # (J/Kg*K) Gas Constant

    Ls = 153, # (deg) Local Siderial Day, 0 = Vernal Equinox, 153 = Viking Lander, 90 = Summer Solstice
    lat = 22.3, # (deg) Latitude

    optical_depth = .5, # Typical on clear days #assumption
    P0 = 669.0, # (pa) Pressure @ Surface Level
    emissGround = .95, #assumption
    albedo = 0.17, #assumption
)

dt = 1. # (s) Time Step for integrating
