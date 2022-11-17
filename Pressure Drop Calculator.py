from Main_Pressure_Drop_Functions import *
import numpy as np
import fluids

class NodeData:
    def __init__(self, deltaP, u, Re, f):
        self.deltaP = deltaP
        self.u = u
        self.Re = Re
        self.f = f

def Pa_to_psi(num):
    return num / 6894.76

tube_network = []
network_data = []

# ------------------------------------------ [Work Here Only] ------------------------------------------ #

# You can set these individually too if you want.
e = 0.000001    # Surface roughness for tubes (meters)
D = 0.0127      # 0.5 inch tube (meters)

# Fluid: Water
rho = 997       # kg/m3
mu = 0.001002   # Dynamic viscosity (Pa s)
m_dot = 0.5     # kg/s
fluid = IncompressibleFluid(rho, mu)

fitting_A = np.pi * D**2 / 4
K = fluids.fittings.K_branch_converging_Crane(D, D, m_dot/rho, m_dot/2/rho)

tube_network.append(RoundTube(e, D, 1))
tube_network.append(Fitting(K, fitting_A))
tube_network.append(RoundTube(e, D, 0.5))
tube_network.append(RoundTube(e, D, 1.5))
tube_network.append(RoundTube(e, D, 3))

# ------------------------------------------------------------------------------------------------------ #

total_deltaP = 0

for item in tube_network:
    u = m_dot / item.A / fluid.rho

    if isinstance(item, Tube):
        deltaP = SP_incomp_deltaP(item, fluid, m_dot)
        reynolds = Re(fluid.rho, u, item.D_h, fluid.mu)
        f = darcy_f(item.e, item.D_h, reynolds)
    else:
        deltaP = item.K * 0.5 * fluid.rho * u**2
        reynolds = None
        f = None

    total_deltaP += deltaP
    network_data.append(NodeData(deltaP, u, reynolds, f))


print(f"Total Pressure Loss: {Pa_to_psi(total_deltaP)} Psi\n")

for i, data in enumerate(network_data):
    if isinstance(tube_network[i], Tube):
        print(f"---------[Tube]---------")
    else:
        print(f"---------[Fitting]---------")

    print(f"Pressure Drop: {Pa_to_psi(data.deltaP)} Psi")
    print(f"Flow Velocity: {data.u} m/s")
    print(f"Reynold\'s Number: {data.Re}")
    print(f"Darcy Friction Factor: {data.f}\n")