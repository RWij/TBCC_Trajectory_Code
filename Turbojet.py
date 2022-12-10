# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 01:31:20 2022

@author: ryant
"""
import math

def Thrust(madot, f, ue, u):
    a = (1+f)*ue
    return madot * (a - u)

# returns To2
def Compr_Inlet_Stag_Temp(mach, Ta):
    a = 1 + (0.5*(1.4-1)*mach*mach)
    return Ta * a

# returns Po2
def Compr_Inlet_Stag_Pres(To2_Ta, Pa):
    nd = 1.0
    y = 1.4
    a = nd * (To2_Ta - 1)
    b = math.pow(1 + a, y / (y - 1))
    return Pa * b

# returns Po3
def Comp_Outlet_Stag_Pres(Prc, Po2):
    return Prc * Po2

# returns To3
def Comp_Outlet_Stag_Temp(Prc, To2):
    nc = 1.0
    y = 1.4
    a = math.pow(Prc, (y - 1) / y)
    b = (1 / nc) * (a - 1)
    return To2 * (1 + b)
    
def fuel_air_ratio(To4_To3, QR, cp, To3):
    a = To4_To3 - 1
    b = (QR / (cp * To3)) - To4_To3
    return a / b

def Turb_Inlet_Pres(Prb, Po3):
    """
    Parameters
    ----------
    Prb : float
        Burner Pressure Ratio (Po4/Po3).
    Po3 : float
        Stagnation Pressure Outside Burner.

    Returns
    -------
    Po4: float
        Stagnation Pressure After Burner

    """
    Po4 = Po3 * Prb
    return Po4

def Turb_Outlet_Stag_Temp(To3, To2, To4):
    """
    Parameters
    ----------
    To3 : float
        Compressor outlet stagnation temperature (K).
    To2 : float
        Compressor inlet stagnation temperature (K).
    To4 : float
        Maximum turbine material temperature (K).

    Returns
    -------
    To5 : float
        Turbine outlet temperature (K).

    """
    
    To5 = To4 - To3 + To2
    return To5

def Turb_Outlet_Stag_Pres(To5_To4, Po4):
    nt = 1.0
    y = 1.4
    a = y / (y - 1)
    b = 1 - ((1.0/nt)*(1 - To5_To4))
    
    print(b)
    if b < 0.0:
        b = 0.0
    Po5 = Po4*math.pow(b, a)
    return Po5

def Exhaust_Velocity(To5, Pa_Po5):
    y = 1.4
    nn = 1.0
    R = 287
    a = y / (y - 1)
    b = 1 - math.pow(Pa_Po5, (y - 1) / y)
    if b < 0.0:
        b = 0.0
    ue = math.pow(2 * nn * a * R * To5 * b, 0.5)
    return ue
    
def TSFC(f, T):
    return f/T

def mass_flow_rate(u, diam, rho):
    area = math.pi*math.pow(diam*0.5,2)
    return area * u * rho

def therm_eff():
    None

def prop_eff(mach, v6_v1, T4_T1, T3_T1):
    y = 1.4
    a = (y - 1)*math.pow(mach, 2)
    b = v6_v1 - 1
    c = (T4_T1) - (T3_T1)
    return a * (b / c)