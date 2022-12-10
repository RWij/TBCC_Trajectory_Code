# -* coding: utf-8 -*-
"""
Created on Tue Apr 12 21:46:34 2022

@author: ryant
"""

import math
import Ramjet
import Turbojet
from atmosphere import atmos

import numpy as np
import openmdao.api as om
#from scipy.optimize import minimize


class JetAnalysisComp(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('num_nodes', types=int)

    def setup(self):
       # inlet conditions
       self.add_input(name='ua', val=0.01, units='m/s', desc='Inlet Velocity')
       self.add_input(name='a', val=340.0, units='m/s', desc='Atmospheric speed of sound')
       self.add_input(name='pa', val=101300, units='N/m**2', desc='Atmospheric Pressure')
       self.add_input(name='Ta', val=288.5, units='K', desc='Atmospheric Temperature')
       self.add_input(name='rho', val=1.225, units='kg/m**3', desc='Atmospheric Density')
       
       # parameters of interest
       self.add_input(name='Prc', val=5.0, desc='Compressor Pressure Ratio')
       self.add_input(name='Prb', val=1.0, desc='Burner Pressure Ratio')
       self.add_input(name='di', val=0.1, units='m', desc='Inlet Diameter')
       self.add_input(name='Fuel', val=0.0, units='kg', desc='Amount of Fuel')
       
       # constant relative to problem
       self.add_input(name='QR', val=42800000.0, units='J/kg', desc='Fuel Heating Ratio')
       self.add_input(name='Tmax', val=2500.0, units='K', desc='Maximum Temperature Limit')
       self.add_input(name='cp', val=1005.0, units='J/(kg*K)', desc='Air Specific Heat at Constant Pressure')

       # output
       self.add_output(name='TSFC', val=0.0, units='kg/(N*s)', desc='Thrust Specific Fuel Consumption')
       self.add_output(name='ue', val=0.0, units='m/s', desc='Outlet velocity')
       self.add_output(name='Thrust', val=0.0, units='N', desc='Thrust')
       self.add_output(name='Fuel_Cons', val=0.0, units='kg', desc='Amount of Fuel Consumed at Instant')
       self.add_output(name='Mach', val=0.0, desc='Inlet/Flight Mach Number')
       
    def compute(self, inputs, outputs):
       ua = inputs['ua']
       a = inputs['a']
       Pa = inputs['pa']
       Ta = inputs['Ta']
       rho = inputs['rho']
       Prc = inputs['Prc']
       QR = inputs['QR']
       To4 = inputs['Tmax']
       cp = inputs['cp']
       Prb = inputs['Prb']
       di = inputs['di']
       Fuel = inputs['Fuel']
       
       mach = ua/a
       
       tsfc = 0.0
       thrust = 0.0
       ue = 0.0 
       f= 0.0
       madot = 0.0
       
       if mach < 1.0:        
           tj_madot = Turbojet.mass_flow_rate(ua, di, rho)
           print(f"inlet mass inflow: {tj_madot} kg/s")
           To2 = Turbojet.Compr_Inlet_Stag_Temp(mach, Ta)
           To2_Ta = To2/Ta
           Po2 = Turbojet.Compr_Inlet_Stag_Pres(To2_Ta, Pa)
           To3 = Turbojet.Comp_Outlet_Stag_Temp(Prc, To2)
           To4_To3 = To4/To3
           tf_f = Turbojet.fuel_air_ratio(To4_To3, QR, cp, To3)
           To5 = Turbojet.Turb_Outlet_Stag_Temp(To3, To2, To4)
           To5_To4 = To5/To4
           Po3 = Turbojet.Comp_Outlet_Stag_Pres(Prc, Po2)
           Po4 = Turbojet.Turb_Inlet_Pres(Prb, Po3)
           Po5 = Turbojet.Turb_Outlet_Stag_Pres(To5_To4, Po4)
           Pa_Po5 = Pa/Po5
           tj_ue = Turbojet.Exhaust_Velocity(To5, Pa_Po5)
           tj_Thrust = Turbojet.Thrust(tj_madot, tf_f, tj_ue, ua)
           tj_tsfc = Turbojet.TSFC(tf_f, tj_Thrust)
               
# =============================================================================
#            print("Turbojet")
#            print(f"ua: {ua}")
#            print(f"Toa: {Ta}  Poa: {Pa}")
#            print(f"To2: {To2}  Po2: {Po2}")
#            print(f"To2/Ta: {To2_Ta}")
#            print(f"To3: {To3}  Po3: {Po3}")
#            print(f"To4: {To4}  Po4: {Po4}")
#            print(f"To4/To3: {To4_To3}")
#            print(f"To5: {To5}  Po5: {Po5}")
#            print(f"To5/To4: {To5_To4}")
#            print(f"tsfc {tj_tsfc}")
# =============================================================================
     
           f = tf_f
           tsfc = tj_tsfc
           thrust = tj_Thrust
           ue = tj_ue
           madot = tj_madot
       
           # assume turbojets work best in the subsonic regime until sonic
           # use turbojet and ramjet to dash towards high mach until
           # though its expected the mach number to switch from tj to rj fully is ~2.2
       elif mach >= 1.0 and mach < 2.2:
            tj_madot = Turbojet.mass_flow_rate(ua, di, rho)
            print(f"inlet mass inflow: {tj_madot} kg/s")
            To2 = Turbojet.Compr_Inlet_Stag_Temp(mach, Ta)
            To2_Ta = To2/Ta
            Po2 = Turbojet.Compr_Inlet_Stag_Pres(To2_Ta, Pa)
            To3 = Turbojet.Comp_Outlet_Stag_Temp(Prc, To2)
            To4_To3 = To4/To3
            tf_f = Turbojet.fuel_air_ratio(To4_To3, QR, cp, To3)
            To5 = Turbojet.Turb_Outlet_Stag_Temp(To3, To2, To4)
            To5_To4 = To5/To4
            Po3 = Turbojet.Comp_Outlet_Stag_Pres(Prc, Po2)
            Po4 = Turbojet.Turb_Inlet_Pres(Prb, Po3)
            Po5 = Turbojet.Turb_Outlet_Stag_Pres(To5_To4, Po4)
            Pa_Po5 = Pa/Po5
            tj_ue = Turbojet.Exhaust_Velocity(To5, Pa_Po5)
            tj_Thrust = Turbojet.Thrust(tj_madot, tf_f, tj_ue, ua)
            tj_tsfc = Turbojet.TSFC(tf_f, tj_Thrust)
            # calculation for pure ramjet so the code knows when to only use
            # ramjet or use tj+rj
            f = tf_f
            tsfc = tj_tsfc
            thrust = tj_Thrust
            ue = tj_ue
            madot = tj_madot
            # if turbojet with afterburner (ramjet), where the exhaust gas from
            # the turbojet stage is reheated
            # compare with propulsive efficiency instead
 
            To4_Toa = To4/To5
            Toa = To5
            rj_f = Ramjet.Fuel_Air_Ratio(To4_Toa, QR, cp, Toa)
            T_madot = Ramjet.Thrust_madot_Ratio(mach, Toa, To4_Toa, rj_f)
            rj_tsfc = Ramjet.TSFC(rj_f, T_madot)
            rj_ue = Ramjet.Exhaust_Velocity(To4, Toa, ua)
            rj_thrust = Ramjet.Thrust(madot, rj_f, rj_ue, ua)
            rj_madot = (1/T_madot) * rj_thrust
            
            print(f"T_madot {T_madot}")
            print(f"rj f {rj_f}")
            print(f"rj tsfc {rj_tsfc}")
             
            tsfc = tsfc + rj_tsfc
            thrust = thrust + rj_thrust
            ue = ue + rj_ue
            f = f + rj_f
            madot = madot + rj_madot
       elif mach >= 2.2:
            To4_Toa = To4/Ta
            Toa = Ta
            rj_f = Ramjet.Fuel_Air_Ratio(To4_Toa, QR, cp, Toa)
            T_madot = Ramjet.Thrust_madot_Ratio(mach, Toa, To4_Toa, rj_f)
            rj_tsfc = Ramjet.TSFC(rj_f, T_madot)
            rj_ue = Ramjet.Exhaust_Velocity(To4, Toa, ua)
            rj_madot = Turbojet.mass_flow_rate(ua, di, rho)
            rj_thrust = Ramjet.Thrust(rj_madot, rj_f, rj_ue, ua)
            madot = rj_madot
            print(f"inlet mass inflow: {madot} kg/s")
            f = rj_f
            tsfc = rj_tsfc
            thrust = rj_thrust
            ue = rj_ue
           
       fuel_consumption = (f * madot)
       if Fuel <= 0.0:
           thrust = 0.0
           fuel_consumption = 0.0
           
       
       print(f"Mach: {mach}, \t Thrust: {thrust/1000} kN, \t ,TSFC: {tsfc}, \t ,ua: {ua} m/s \t ue: {ue} m/s \t ,f: {f}") 
                      
       outputs['TSFC'] = tsfc
       outputs['ue'] = ue
       outputs['Thrust']= thrust
       outputs['Fuel_Cons'] = fuel_consumption
       outputs['Mach'] = mach
           
    
class TrajectoryEOM2D(om.ExplicitComponent):
    """
    Computes the position and velocity equations of motion using a 2D flight path.
    """
    def initialize(self):
        self.options.declare('num_nodes', types=int)
                
    def setup(self):
        self.add_input(name='m', val=1.0, units='kg',
                       desc='total Mass')
        self.add_input(name='v', val=0.01, units='m/s',
                       desc='Velocity magnitude')
        self.add_input(name='T', val=0.0, units='N',
                       desc='thrust')
        self.add_input(name='alpha', val=0.0, units='rad',
                       desc='angle of attack')
        self.add_input(name='dt', val=0.2, units='s',
                       desc='Time Increments')

        self.add_output(name='v_dot', val=0.0, units='m/s**2',
                        desc='rate of change of velocity magnitude')
        self.add_output(name='hchng_dot', val=0.0, units='m/s',
                        desc='rate of change of vertical velocity')
        self.add_output(name='rchng_dot', val=0.0, units='m/s',
                        desc='rate of change of horizontal velocity')

    def setup_partials(self):
        self.declare_partials(['hchng_dot', 'rchng_dot'], ['T', 'm','alpha'])

    def compute(self, inputs, outputs):
        g =  9.80665 # m/s^2
        m = inputs['m']
        dt = inputs['dt']
        T = inputs['T']
        alpha = inputs['alpha']
        
        h_2dot = ((T * np.sin(alpha)) / m) - g
        r_2dot = (T * np.cos(alpha)) / m
        hchng_dot = h_2dot * dt
        rchng_dot = r_2dot * dt
        outputs['v_dot'] = math.pow(math.pow(h_2dot,2) + math.pow(r_2dot,2), 0.5)
        outputs['hchng_dot'] = hchng_dot
        outputs['rchng_dot'] = rchng_dot
        
        
class TbccODE(om.Group):

    def setup(self):
        self.add_subsystem(name='jetanal',
                           subsys=JetAnalysisComp(),
                           promotes=['*'])

        self.add_subsystem(name='eom',
                           subsys=TrajectoryEOM2D(),
                           promotes=['*'])

        self.connect('Thrust', 'T')
        

def eval_range(design_init, prob, complex_step=False):
     h = 0.0
     r = 0.0
     hdot = 0.0
     rdot = 0.0
     dt = 0.1
     t = 0.0
     [_,_,_,_,a] = atmos(h)
     Mach = 0.2 # typical takeoff Mach
     ua = a * Mach
     emass = 1600  # empty engine mass (tj and ramjet)
     fuel_mass = 2000.0 # fuel mass at takeoff
     in_mission = True # only false when it crashes from lack of fuel after takeoff during glide descent

     if complex_step:
        prob.set_complex_step_mode(True)

     while in_mission:
        print("=====================================")
        print(f"Time {t} sec.     Mass {fuel_mass+emass} kg")
        [T, P, rho, _, a] = atmos(h)

        # Set values
        prob.set_val('ua', ua, units='m/s')
        prob.set_val('v', ua, units='m/s')
        prob.set_val('a', a, units='m/s')
        prob.set_val('pa', P, units='N/m**2')
        prob.set_val('Ta', T, units='K')
        prob.set_val('rho', rho, units='kg/m**3')
        prob.set_val('m', fuel_mass + emass, units='kg')     
        prob.set_val('Fuel', fuel_mass, units='kg')
        
        prob.set_val(name='Prc', val=design_init)
        prob.set_val(name='Prb', val=12.0)
        prob.set_val(name='di', val=1.0, units='m')
        
        Mach = ua/a
        
        # Run the model
        prob.run_model()

        # Extract rates
        v_dot = prob.get_val('v_dot')
        hchng_dot = prob.get_val('hchng_dot')
        rchng_dot = prob.get_val('rchng_dot')
        f_cons = prob.get_val('Fuel_Cons')
        Mach = prob.get_val('Mach')
                    
        h_last = h
        r_last = r
            
        # Euler Integration
        fuel_mass = fuel_mass - (dt * f_cons)
        ua = ua + (dt * v_dot)
        hdot = hdot + hchng_dot
        rdot = rdot + rchng_dot
        h = h + (dt * hdot)
        r = r + (dt * rdot)
        t += dt
        
        if Mach >= 4.0:
            ua = 4.0*a
            
        if h >= 20000:
            prob.set_val('alpha', 1.0*(math.pi/180.0), units='rad')
            h = 20000
        
        print(f"Height: {h/1000.0} km\t Range: {r/1000.0} km \t fuel consumed: {f_cons*dt} kg")
        print("=====================================") 
        if h <= -0.01:
            in_mission = False
    
     # Linear interpolation between last two points to get the landing point accurate.
     r_final = r_last + (r - r_last) * h_last / (h_last - h)
    
     if complex_step:
         prob.set_complex_step_mode(False)
    
     return -r_final


def gradient_range(init, prob):
    """
    Uses complex step to compute gradient of range wrt initial angle.

    Parameters
    ----------
    init : float
        Initial conditions
    prob : <Problem>
        OpenMDAO problem that contains the equations of motion.

    Returns
    -------
    float
        Derivative of range wrt pressure ratio
    """
    step = 0.1
    dr_dgam = eval_range(init + step, prob, complex_step=False)
    return dr_dgam.imag / step

       
if __name__ == "__main__":   
    prob = om.Problem(model=TbccODE())
    prob.setup(force_alloc_complex=False)
    
    ###### Constants ######
    prob.set_val('QR', 42800.0, units='kJ/kg')
    prob.set_val('Tmax', 2500.0, units='K')
    prob.set_val('cp', 1005.0, units='J/(kg*K)')    
    prob.set_val('alpha', 75.0*(math.pi/180), units='rad')
    
    
    #Prc = np.linspace(4, 40, 20)
    #Prb = np.linspace(4, 40, 20)
    #Diam_inlet = np.linspace(0.05, 0.5, 20) # meters
    ####################
  #  to_optimize = False
    
#    if to_optimize:
        #result = minimize(eval_range, 30.0,
        #               method='SLSQP',
    #                   jac=gradient_range,
    #                   args=(prob))
     #   print(result['x'])
    #else:
        #result = eval_range(45.0, prob)
result = eval_range(13.0, prob)

