exec(open("/home/zihanw8/scratch/dynamics_2020_07/ModulesWithLags.py", encoding = 'utf-8').read())
exec(open("/home/zihanw8/scratch/dynamics_2020_07/parameters.py", encoding = 'utf-8').read())
import numpy as np
from math import *
import random
import itertools
from scipy.optimize import root
from scipy import stats
from scipy import special
import copy
import time
import collections
import pickle
import sys

Nr = int(sys.argv[1])
# index of simulation is marked as sys.argv[2]
lagTime = 0.1
Nb = factorial(Nr)*Size
yields_list = 0.5*np.ones([Nb, Nr])
growth_rate_list  = []
for i in range(Nb):
    growth_rate_list.append([max(np.random.normal(g_mu, g_sigma), theta) for i in range(Nr)])
invlist = list(range(Nb))
random.shuffle(invlist)

inputFileName = '/home/zihanw8/scratch/dynamics_2020_07/input_'+'Nr='+sys.argv[1]+'_index='+sys.argv[2]+'_balanced.pkl'
pickle.dump({'Nr': Nr, 'G': growth_rate_list, 'invlist': invlist, 'Y': yields_list}, open(inputFileName, "wb" ))

details = {'res_begin':[], 't_info':[], 'res_left':[], 'bug_info':[], 'round_idx':[0]}
Res = [1.0 for i in range(Nr)] # for balanced nutrients
# Res = [Nr] + [Nr/100 for i in range(Nr-1)] # for unbalanced nutrients
preference_list = list(itertools.permutations(range(Nr), Nr))
start=time.time()
movestep=0
details = {'res_begin':[], 't_info':[], 'res_left':[], 'bug_info':[], 'round_idx':[0]}
system = {'res_available': np.heaviside(Res, 0), 'res_concentration': [i for i in Res], 'bug_available': [0 for i in range(Nb)], 'bug_density': [0 for i in range(Nb)]}
# ext_list: list of bugs not in the community
old_ext_list = invlist
#print(old_ext_list)
system, new_ext_list = round_robin_invade(system, old_ext_list, growth_rate_list)
new_ext_list = [i for i in invlist if i in new_ext_list]
count = 0
cycled = False
while new_ext_list != old_ext_list and count < 10:
    old_ext_list = new_ext_list
    survivors, concent, pref_list, growth = output(system)
    #print(survivors)
    system, new_ext_list = round_robin_invade(system, old_ext_list, growth_rate_list)
    count = count + 1
    if count > 6:
        cycled = True
        print('Could have been cycling')
#print(survivors)
#print(output(system))
details['round_idx'] = details['round_idx'][1:]
end = time.time()
print(end-start)
#print(getMarketHistory(details))

outputFileName = '/home/zihanw8/scratch/dynamics_2020_07/outputs_'+'Nr='+sys.argv[1]+'_index='+sys.argv[2]+'_balanced.pkl'
pickle.dump({'details': getMarketHistory(details), 'output': output(system), 'cycled' : cycled}, open(outputFileName, "wb" ))
