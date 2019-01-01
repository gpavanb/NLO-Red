# Storage directory
directory = 'FLAMES'

# Reference palette
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
enable_dc = False
tolerance_dc = 0.0125

# Ai cases
# Pressure
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
x0 = [0.02244707,0.02122669, 0.02100455, 0.0215207, 0.0215207, 0.02764953,
 0.02509764, 0.02509764, 0.02665645, 0.0350278, 0.0339066, 0.0339066,
 0.03000983, 0.03042798, 0.03962287, 0.04125298, 0.04125025, 0.04858254,
 0.04732721, 0.04619811, 0.04351466, 0.05390394, 0.060274, 0.06022112,
 0.07117585, 0.07117671]
