exec(open("Modules.py", encoding = 'utf-8').read())
exec(open("parameters.py", encoding = 'utf-8').read())
import numpy as np
from math import *
import random
import itertools
from scipy.optimize import root
from scipy import stats
from scipy import special
import copy
import time
# import matplotlib.pyplot as plt
# from matplotlib.patches import Circle
# import matplotlib.patches as mpatches
# import seaborn as sns
import collections
import pickle
import sys
import multiprocessing

def model(Nr_in, indices):

	global Nr
	global Nb
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

	Nr = Nr_in
	Nb = factorial(Nr)*Size
	yields_list = 0.5*np.ones([Nb, Nr])
	growth_rate_list  = []
	for i in range(Nb):
	    growth_rate_list.append([max(np.random.normal(g_mu, g_sigma), theta) for i in range(Nr)])
	invlist = list(range(Nb))
	random.shuffle(invlist)

	inputFileName = 'input_'+'Nr='+str(Nr)+'_index='+str(indices)+'_balanced.pkl'
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
	    if count > Nr:
	        cycled = True
	        print('This looks cycled')
	details['round_idx'] = details['round_idx'][1:]
	end = time.time()

	outputFileName = 'outputs_'+'Nr='+str(Nr)+'_index='+str(indices)+'_balanced.pkl'
	pickle.dump({'details': getMarketHistory(details), 'output': output(system), 'cycled' : cycled}, open(outputFileName, "wb" ))

	return end-start


if __name__ == '__main__':
	N_indices = 20
	cores = multiprocessing.cpu_count()
	pool = multiprocessing.Pool(processes=1)
	x1 = [Nr for i in range(N_indices)]
	x2 = list(range(N_indices))
	tasks = [(Nr_in, indices) for Nr_in in x1 for indices in x2]
	print(pool.starmap(model,tasks))


