"""Physical constants (SI, CODATA 2018) and the practical units used throughout.

Practical units, matching how solar wind data are published:
    magnetic field  nT        number density  cm^-3
    velocity        km/s      temperature     K
"""

MU0 = 1.25663706212e-6  # vacuum permeability, N A^-2
KB = 1.380649e-23  # Boltzmann constant, J K^-1
MP = 1.67262192369e-27  # proton mass, kg

NT = 1e-9  # nT -> T
CC = 1e6  # cm^-3 -> m^-3
KMS = 1e3  # km/s -> m/s

GAMMA = 5.0 / 3.0  # adiabatic index of a monatomic ideal gas
