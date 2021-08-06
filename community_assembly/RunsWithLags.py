exec(open("ModulesWithLags.py", encoding = 'utf-8').read())
exec(open("Parameters.py", encoding = 'utf-8').read())
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
import multiprocessing
from multiprocessing import Pool, freeze_support

def model(Nr_in, indices):

    global Size
    global D
    global yields_list
    global growth_rate_list
    global invlist
    global details
    global Res
    global preference_list
    global system
    global movestep

    global g_mu
    global g_var
    global scale
    global theta_low
    global theta_high_more
    global lagTime
    global Nb
    global Nr

    Nr=Nr_in
    Nb = factorial(Nr)*Size
    yields_list = 0.5*np.ones([Nb, Nr])
    growth_rate_list  = []
    for i in range(Nb):
        growth_rate_list.append([max(np.random.normal(g_mu, g_sigma), theta) for i in range(Nr)])
    invlist = list(range(Nb))
    random.shuffle(invlist)

    inputFileName = 'input_'+str(indices)+'_lag='+str(lagTime)+'.pkl'
    pickle.dump({'Nr': Nr, 'G': growth_rate_list, 'invlist': invlist, 'Y': yields_list}, open(inputFileName, "wb" ))

    details = {'res_begin':[], 't_info':[], 'res_left':[], 'bug_info':[], 'round_idx':[0]}
    Res = [1.0 for i in range(Nr)] # for balanced nutrients
    # Res = [Nr] + [Nr/100 for i in range(Nr-1)] # for unbalanced nutrients
    Res = [i/sum(Res)*Nr for i in Res]
    preference_list = list(itertools.permutations(range(Nr), Nr))
    start=time.time()
    movestep=0
    details = {'res_begin':[], 't_info':[], 'res_left':[], 'bug_info':[], 'round_idx':[0]}
    system = {'res_available': np.heaviside(Res, 0), 'res_concentration': [i for i in Res], 'bug_available': [0 for i in range(Nb)], 'bug_density': [0 for i in range(Nb)]}
    # ext_list: list of bugs not in the community
    old_ext_list = invlist
    system, new_ext_list = round_robin_invade(system, old_ext_list, growth_rate_list)
    new_ext_list = [i for i in invlist if i in new_ext_list]
    count = 0
    cycled = False
    while new_ext_list != old_ext_list and count < 10:
        old_ext_list = new_ext_list
        survivors, concent, pref_list, growth = output(system)
        system, new_ext_list = round_robin_invade(system, old_ext_list, growth_rate_list)
        count = count + 1
        if count > 8:
            cycled = True
            print('A cycle might have occured')
    details['round_idx'] = details['round_idx'][1:]
    end = time.time()
    outputFileName = 'outputs_'+str(indices)+'_lag='+str(lagTime)+'.pkl'
    pickle.dump({'details': getMarketHistory(details), 'output': output(system), 'cycled' : cycled}, open(outputFileName, "wb" ))
    return (end-start)

if __name__ == '__main__':
	freeze_support()
	cores = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(processes=9)
	x1 = [Nr for i in range(N_indices)]
	x2 = list(range(N_indices))
	tasks = [(Nr_in, indices) for Nr_in, indices in zip(x1,x2)]
	pool.starmap(model,tasks)