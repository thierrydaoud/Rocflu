"""Scales between the dimensionless problems and Rocflu's SI input: 1e5 Pa, 1 kg/m3, 1 m."""
import math

P_SCALE = 1.0e5
RHO_SCALE = 1.0
U_SCALE = math.sqrt(P_SCALE / RHO_SCALE)  # m/s
T_SCALE = 1.0 / U_SCALE                   # s, for a 1 m domain
CP = 1004.64                              # J/(kg K); with gamma fixed it only sets the temperature


def si(state):
    """(rho, u, p) dimensionless -> SI."""
    return state[0] * RHO_SCALE, state[1] * U_SCALE, state[2] * P_SCALE
