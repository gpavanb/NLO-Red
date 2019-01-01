# Storage directory
directory = 'FLAMES'

# Reference palette
# Stanford A
#palette = ['MCYC6','C7H8','C6H6','IC8H18','NC12H26']
#test_comp = [0.1, 0.1, 0.01, 0.055, 0.735]
# Dryer
test_comp = [0.404,0.295,0.073,0.228]
palette = ['NPBENZ','NC12H26','IC8H18','TMBENZ']

# Species list
spec_list = 'data/POLIMI_List.txt'

# Cantera mechanism
mechanism = 'mech/POLIMI_TOT.xml'

# Species to exclude (NOx chemistry for example)
species_exclude_init = ()

# Major species, not in the optimization loop
species_major = ('N2', 'AR', 'CH3OCH3', 'CO2', 'H2O',
                 'O2', 'H2', 'CO', 'OH', 'HO2')

# Fuel and N2-O2 ratio
fuel = 'CH3OCH3'
n2_o2_ratio = 3.76

# Threshold value for Zero
threshold = 0.05

# Error tolerances
tolerance_mw = 0.05
tolerance_hc = 0.005

# Sum relaxation
sum_relax = 1e-3

# Active subspace
use_active_subspace = True

# Distillation curve cases
enable_dc = True
tolerance_dc = 0.0125

# Ai cases
# Pressure
enable_idt = False
P_ai = [1e5]
# Temperature
# T_ai = [1200.0, 1400.0, 1600.0]
T_ai = [1200.0]
# Equivalence ratio
phi_ai = [1.0]
# phi_ai = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
#          1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
# Error tolerance
tolerance_ai = 0.20

# Flame cases
P_fl = []
T_fl = []
phi_fl = []

# Flame cases
# Pressure
# P_fl = [1e5]
# # Temperature
# T_fl = [300.0]
# # Equivalence ratio
# phi_fl = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
#           1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
# Error tolerance
tolerance_fl = 0.05

# Initial beta vector
# x0 = 
