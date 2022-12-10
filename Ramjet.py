# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 01:20:31 2022

Functions for an Ideal Ramjet 

@author: ryant
"""
import math

def Thrust(madot, f, ue, u):
    return madot*((1+f)*ue - u)

def Isentropic_Total_Pressure(mach, pres):
    a = 1.4 / (1.4 - 1)
    b = (1.4 - 1) * 0.5
    c = 1 + (b*mach*mach)
    return pres * math.pow(c, a)

def Exhaust_Velocity(To4, Toa, ua):
    return ua * math.pow(To4/Toa,0.5)

def Fuel_Air_Ratio(To4_Toa, QR, cp, Toa):
    a = To4_Toa - 1
    b = (QR/(cp*Toa)) - To4_Toa
    return a/b

def Thrust_madot_Ratio(mach, Toa, To4_Toa, f):
    a = mach * math.pow(1.4*287*Toa, 0.5)
    b = (1+f)*math.pow(To4_Toa, 0.5)
    c = 1 + ((1.4 - 1) * 0.5*mach*mach)
    d = (b * math.pow(c, -0.5)) - 1
    return a * d

def TSFC(f, T_madot):
    return f/T_madot
    
def mfdot(TSFC, Thrust):
    return TSFC*Thrust

def mass_flow_rate(u, diam, rho):
    area = math.pi*math.pow(diam*0.5,2)
    return area * u * rho

def therm_eff(prc):
    return 1 - (1/prc)

def prop_eff():
    return 