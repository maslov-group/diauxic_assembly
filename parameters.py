Nr = 4

# distribution of growth rates, with symmetric nutrients
theta = 1e-1 # hard minimum of growth rates
g_mu  = 0.25
g_sigma= 0.05
Size = 480 # 480 different species for each preference order
lagTime = 0
# dynamic parameters
D = 0.01 # (inverse) dilution factor
T_dilute = 24 # time interval between dilutions
dilute_to_steady = 10 # #(dilution) between invasions
b0 = 1e-6 # density of bug when introduced/initial
b_threshold = 1e-6 # extinction density
c_threshold = 1e-9 # concentration threshold