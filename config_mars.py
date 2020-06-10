balloon_properties = dict(
    shape = 'sphere',
    d = 20.0,               # (m) Diameter of Sphere Balloon
    mp = 10.,                # (kg) Mass of Payload

    cp = 900.0,             #(J/(kg K)) Specific heat of envelope material Germanium : 320, Aluminum: 900
    emissEnv = .03,         # Emisivity of Enevelope
    absEnv = .6,            # Absorbiviy of envelope
)

control_properties = dict(
    vent = 0.0,             # (kg/s) Vent Mass Flow Rate
    alt_sp = 2000.,         # (m) Altitude Setpoint
    v_sp = 0.,              # (m) Altitude Setpoint, Not Implemented right now
)

mars_properties = dict(
    Cp_co2 = 735.0,         # (J/Kg*K)  Specifc Heat Capacity, Constant Pressure
    Cv_co2 = 657.0,         # (J/Kg*K)  Specifc Heat Capacity, Constant Volume
    Rsp_co2 = 188.92,       # (J/Kg*K) Gas Constant

    Ls = 90.,                # (deg) Local Siderial Day, 0 = Vernal Equinox, 153 = Viking Lander, 90 = Summer Solstice
    lat = 22.3,             # (deg) Latitude

    optical_depth = .5,     # Typical on clear days #assumption
    P0 = 669.0,             # (Pa) Pressure @ Surface Level
    emissGround = .95,      # Emissivity of Ground
    albedo = 0.17,          # Albedo of Ground
)

dt = 1. #Do not change the stepsize for now, Only implemented for 1s.
