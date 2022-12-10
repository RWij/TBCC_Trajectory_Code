# -*- coding: utf-8 -*-
"""
Modified on Wed Apr 13 02:10:15 2022
By Ryan Wijaya 
for use in metric
"""

import math

def atmos(h):
    """
    1976 US Standard Atmosphere
    http://www.atmosculator.com/The#20Standard#20Atmosphere.html?
    David Pate and Michael Patterson

    Returns
        A list with the following:
          
    Repurposed from AE 4802: Configuration Aerodynamics and Flight Performance
    by Brian German

    Parameters
    ----------
    h : float
        altitude in feet.

    Returns
    -------
    temp:  float
        temperature in K
    pres: float
        pressure in Pa
    rho: float
        density in kg/m^3
    mu: float
        viscosity in kg/(m-s)
    a: float
        sound speed m/s assuming R=1716 ft-lb/(slugs-°R) and gamma = 1.4
    """
    
    R = 1716;
    gamma = 1.4;    
    T_SL = 518.69;      #°R
    T0 = 491.6;         #°R
    p_SL = 2116.2;      #lb/ft^2
    rho_SL = .0023769;  #slug/ft^3
    mu0 = 3.58394051e-7;    #slug/(ft-s)
    # 3.7373e-7;    #slug/(ft-s)
    S = 199;            #Sutherland's Constant °R
    theta = 0
    delta = 0
    sigma = 0
    
    h = h * 3.281
    
    if h <= 36809:
        theta = 1 - h / 145442;
        delta = math.pow(1 - h / 145442, 5.255876);
        sigma = math.pow(1 - h / 145442, 4.255876);
        
    
    # Isothermal
    elif h > 36089 and h <= 65617:
        theta = 0.751865;
        delta = 0.223361 * math.exp(-(h-36089)/20806);
        sigma = 0.297076 * math.exp(-(h-36089)/20806);
        
    
    # (Inversion)
    elif h > 65617 and h <= 104987:
        theta = 0.682457 + h / 945374;
        delta = math.pow(0.988626 + h / 652600, -34.16320);
        sigma = math.pow(0.978261 + h / 659515, -35.16320);
    
    
    # (Inversion)
    elif h > 104987 and h <= 154199:
        theta = 0.482561 + h / 337634;
        delta = math.pow(0.898309 + h / 181373, -12.20114);
        sigma = math.pow(0.857003 + h / 190115, -13.20114);
    
    
    # (Isothermal)
    elif h > 154199 and h <= 167323:
        theta = 0.939268;
        delta = 0.00109456 * math.exp(-(h-154199)/25992);
        sigma = 0.00116533 * math.exp(-(h-154199)/25992);
    
    
    elif h > 167323 and h <= 232940:
        theta = 1.434843 - h / 337634;
        delta = math.pow(0.838263 - h / 577922, 12.20114);
        sigma = math.pow(0.798990 - h / 606330, 11.20114);
    
    
    elif h > 232940 and h <= 278386:
        theta = 1.237723 - h / 472687;
        delta = math.pow(0.917131 - h / 637919, 17.08160);
        sigma = math.pow(0.900194 - h / 649922, 16.08160);
    
    
    temp = theta * T_SL;
    pres = delta * p_SL;
    rho = sigma * rho_SL;
    mu = mu0 * math.pow(temp/T0,3/2) * (T0 + S) / (temp + S);
    a = math.pow(gamma * R * temp,0.5);
    
    # convert to metric
    temp = temp * 5/9
    pres = (pres / 20.88543)*1000
    rho = rho * 14.594/(0.3048*0.3048*0.3048)
    mu = mu * 14.594 * 3.281
    a = a / 3.281
        
    return [temp, pres, rho, mu, a]
