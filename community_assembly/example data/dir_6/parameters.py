Nr=4
Size = 480

# distribution of growth rates, with symmetric nutrients
theta = 1e-1 # hard minimum of growth rates
g_mu  = 0.5
g_sigma= 0.1

# dynamic parameters
D = 0.01 # (inverse) dilution ratio
T_dilute = 12 # time interval between dilutions
dilute_to_steady = 10 # #(dilution) between invasions
b0 = 1e-6 # density of bug when introduced/initial
b_threshold = 1e-6 # extinction density
c_threshold = 1e-9 # concentration threshold

