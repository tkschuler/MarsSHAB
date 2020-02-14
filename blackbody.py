import math
from scipy import optimize

SB = 5.670373E-8


r = 2

alpha = .6
emmis = .03
a_solar = math.pi*pow(r,2)
a = 4*math.pi*pow(r,2)
SF = 500


Pin = alpha*a_solar*SF


def f(T_s):
    return Pin-emmis*SB*a*(pow(T_s,4))

T_s = optimize.newton(f, 400)
print T_s
#print Pin-emmis*SB*a*(pow(T_s,4))
