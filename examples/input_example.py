# Storage directory
directory = 'FLAMES'

# Reference palette
#palette = ['MCYC6','C7H8','C6H6','IC8H18','NC12H26']
# Dryer
palette = ['NPBENZ','NC12H26','IC8H18','TMBENZ']

# Reference composition
#test_comp = [0.1, 0.1, 0.01, 0.055, 0.725]
test_comp = [0.404,0.295,0.073,0.218]

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

# Distillation curve cases
enable_dc = True
tolerance_dc = 0.0125

# Ai cases
# Pressure
P_ai = []
# Temperature
# T_ai = [1200.0, 1400.0, 1600.0]
T_ai = []
# Equivalence ratio
phi_ai = []
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
# x0 = [9.99999944e-01, 9.99999878e-01, 9.99999904e-01, 9.99999960e-01
# , 5.37093170e-01, 5.08976098e-01, 6.01165069e-08, 6.86705244e-08
# , 6.47897171e-08, 6.17028363e-02, 9.99999914e-01, 3.20745868e-08
# , 9.99999951e-01, 1.02802446e-03, 1.36235472e-07, 9.99999694e-01
# , 8.76082885e-08, 2.75813114e-02, 5.54320308e-08, 1.07509611e-07
# , 7.68159552e-08, 9.66355456e-08, 5.96294581e-08, 1.07365982e-07
# , 8.39643376e-08, 9.60614566e-08, 1.03017202e-07, 9.21567832e-08
# , 9.99999793e-01, 7.08789395e-08, 1.06592199e-07, 4.66776158e-08
# , 9.47730344e-08, 1.47151161e-07, 1.10629874e-07, 1.07826961e-07
# , 7.89858822e-08, 1.05792538e-07, 1.01698189e-07, 9.23493802e-08
# , 9.84981048e-08, 1.01632182e-07, 8.47171651e-08, 1.06656556e-07
# , 1.05745840e-07]
