# Author: Jiyansh Agarwal

import numpy as _np
from scipy.optimize import fsolve as _fsolve

# Notes:
# Straight Pipe = SP
# All pipes are assumed to be constant cross sectional area
# All units must be SI

# ------------------------ [Classes] ------------------------ #

# e = pipe roughness factor (NOT relative roughness)
# D_h = pipe hydraulic diameter
# A = pipe cross-sectional area
# L = pipe length
class Tube:
    def __init__(self, e, D_h, A, L):
        self.e = e
        self.D_h = D_h
        self.A = A
        self.L = L

# e = pipe roughness factor (NOT relative roughness)
# D = pipe diameter
# L = pipe length
class RoundTube(Tube):
    def __init__(self, e, D, L):
        super().__init__(e, D, _np.pi * D**2 / 4, L)


class Fitting:
    def __init__(self, K, A):
        self.K = K
        self.A = A

# rho = fluid density
# mu = dynamic viscosity
# g = ratio of specific heats, gamma
# Rs = Specific gas constant
class CompressibleFluid:
    def __init__(self, rho, mu, g, Rs):
        self.rho = rho
        self.mu = mu
        self.g = g
        self.Rs = Rs

# rho = fluid density
# mu = dynamic viscosity
class IncompressibleFluid:
    def __init__(self, rho, mu):
        self.rho = rho
        self.mu = mu

# ------------------------ [Functions] ------------------------ #

# rho = density
# u = velocity
# D_h = pipe hydraulic diameter
# mu = dynamic viscosity
def Re(rho, u, D_h, mu):
    return rho * u * D_h / mu

# e = pipe roughness factor (NOT relative roughness)
# D_h = pipe hydraulic diameter
# Reynold's number of the flow
def darcy_f(e, D_h, Re):
    if Re < 2100:       # Laminar Flow
        return 64 / Re
    else:               # Turbulent Flow if Re > 4000
        if Re < 4000:
            print("Warning: Transition flow, friction factor may not be completely accurate.")
        
        # Solving using Colebrook-White Equation
        #func = lambda f : 1 / _np.sqrt(f) + 2 * _np.log10(e / 3.7 / D_h + 2.51 / Re / _np.sqrt(f))
        #return _fsolve(func, 0.1)[0]

        # Solving using Swamee–Jain Equation (Full-flow circular tube)
        return 0.25 / (_np.log10(e / 3.7 / D_h + 5.74 / Re**0.9) ** 2)


# m_dot = mass flow rate
def SP_incomp_deltaP(tube: Tube, incomp_fluid: IncompressibleFluid, m_dot):
    u = m_dot / incomp_fluid.rho / tube.A
    reynolds_num = Re(incomp_fluid.rho, u, tube.D_h, incomp_fluid.mu)
    f = darcy_f(tube.e, tube.D_h, reynolds_num)
    
    # Darcy–Weisbach equation
    return (f * tube.L / tube.D_h) * 0.5 * incomp_fluid.rho * u**2

# ------------------------ [Compressible Flow] ------------------------ #

# Everywhere g = ratio of specific heats, gamma

# Fanno flow equation. M is the entrance mach number.
def fl_d_max(M, g):
    return ((g+1)/2/g) * _np.log(M**2 * (g+1) / (2 + (g-1) * M**2)) - (1 - 1 / M**2) / g

# Ratio of P1 to P2 for any two points in duct. M1 and M2 are the corresponding mach numbers.
def fanno_P1_P2(M1, M2, g):
    return (M2 / M1) * _np.sqrt((2 + (g-1) * M2 ** 2) / (2 + (g-1) * M1 ** 2))

# Ratio of T1 to T2 for any two points in duct. M1 and M2 are the corresponding mach numbers.
def fanno_T1_T2(M1, M2, g):
    return (2 + (g-1) * M2 ** 2) / (2 + (g-1) * M1 ** 2)

# Ratio duct area to critical duct area (location where M=1) for isentropic flow. M is the mach number.
def isentropic_A_A_star(M, g):
    return ((2 / (g+1)) * (1 + 0.5 * (g-1) * M**2)) ** ((g+1)/2/(g-1)) / M



# This applies fanno flow (adiabatic, constant-area duct with friction) to ideal gases only. 
# The given conditions are at the entrance of the tube and the function returns the conditions at the pipe exit in a list as follows:
#   [ Exit Temperature, Exit Pressure, Exit Mach number]

# m_dot = mass flow rate
# T1 = inlet temperature
# P1 = inlet pressure

def SP_Fanno_Flow(tube: Tube, comp_fluid: CompressibleFluid, m_dot, T1, P1):
    u1 = m_dot / comp_fluid.rho / tube.A
    a1 = _np.sqrt(comp_fluid.g * comp_fluid.Rs * T1)    # Speed of sound
    M1 = u1 / a1

    reynolds_num = Re(comp_fluid.rho, u1, tube.D_h, comp_fluid.mu)
    f = darcy_f(tube.e, tube.D_h, reynolds_num)

    fl_d_2 = fl_d_max(M1, comp_fluid.g) - f * tube.L / tube.D_h  # This indicates the "remaining" distance till M=1 is reached

    # Now numerically find the mach number that corresponds to fl_d_2
    func = lambda M2 : fl_d_max(M2, comp_fluid.g) - fl_d_2
    M2 = _fsolve(func, M1)[0]

    P2 = P1 / fanno_P1_P2(M1, M2, comp_fluid.g)
    T2 = T1 / fanno_T1_T2(M1, M2, comp_fluid.g)

    return [T2, P2, M2]