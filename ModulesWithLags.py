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

# functions involved in the serial dilution simulation
def dilute_check(system, growth_rate_list, lagtime = 0, ignore_this_run = False, ReturnInvasionChance = True):
    
    # t_points is the end of each nutrient's depletion, as well as 0 and T_dilution(period).
    t_points = []
    # c_points is the bug densities of each t point.
    c_points = []
    # r_points is the available nutrients at each t point.
    r_points = []
    #conc_points is the concentration of each nutrient.
    conc_points = []
    # inv_points is the measure of conflict or invasibility at each t point. 
    inv_points = []
    t_switch = 0
    t_points.append(t_switch)
    survivors, concent, pref_list, growth = output(system)
    c_points.append(concent)
    conc_points.append([i for i in system['res_concentration']])
    lag_flag = [0 for i in range(Nb)]
    lag_begins = [0 for i in range(Nb)]
    #list of resource in use for each consumer
    use = [0 for i in range(Nb)]
    while t_switch < T_dilute:
        a = [i for i in system['res_available']]
        r_points.append(a)
        # list of consumers for each resource
        consumer = [[] for i in range(Nr)]
        for i in range(Nb):
            # if all are depleted, bug still uses its least preferred nutrient
            while system['res_available'][preference_list[i//Size][use[i]]] < 1 and use[i] < Nr - 1: # when there's no current resource left
                use[i] = use[i] + 1 # goto next resource
                lag_flag[i] = 1 # start the lag
                lag_begins[i] = t_switch # and mark the starting point. Notice: lags do superposition with each other.
            if system['bug_available'][i] > 0:
                consumer[preference_list[i//Size][use[i]]].append(i)
        # find the earliest depleted resource
        t_dep = T_dilute - t_switch
        for i in range(Nr):
            if system['res_available'][i] > 0:
                def remain(t):
                    # S = c - sum(B0(e^gt-1)/Y)
                    return system['res_concentration'][i] - sum([system['bug_density'][j]*(exp(t*growth_rate_list[j][i]) - 1)/yields_list[j][i] for j in consumer[i] if lag_flag[j] == 0])
                # a little bit tricky here...
                t_i = root(remain, T_dilute + 1).x[0]
                if t_i < t_dep:
                    t_dep = t_i
        for i in range(Nb):
            if lag_flag[i] == 1:
                t_i = lag_begins[i] + float(lagTime) - t_switch
                if t_i < t_dep:
                    t_dep = t_i
        # update the system according to this t_dep

        ResInUse = sum([int(i!=[]) for i in consumer])
        ResAvailable = sum(system['res_available'])
        if(ResAvailable!=0):
            inv_points.append(1-ResInUse/ResAvailable)
        t_switch = t_switch + t_dep
        temp_bug_density = [i for i in system['bug_density']]
        for i in range(Nb):
            if system['res_available'][preference_list[i//Size][use[i]]] > 0 and lag_flag[i] == 0:
                temp_bug_density[i] = system['bug_density'][i]*exp(growth_rate_list[i][preference_list[i//Size][use[i]]]*t_dep)
        for i in range(Nr):
            if system['res_available'][i] > 0:
                system['res_concentration'][i] = system['res_concentration'][i] - sum([system['bug_density'][j]*(exp(t_dep*growth_rate_list[j][i]) - 1)/yields_list[j][i] for j in consumer[i] if lag_flag[j] == 0])
                if system['res_concentration'][i] < c_threshold:
                    system['res_available'][i] = 0
        for i in range(Nb):
            system['bug_density'][i] = temp_bug_density[i]
        t_points.append(t_switch)
        survivors, concent, pref_list, growth = output(system)
        c_points.append(concent)
        conc_points.append([i for i in system['res_concentration']])
        for i in range(Nb):
            if lag_flag[i] == 1 and lag_begins[i] + float(lagTime) == t_dep + t_switch:
                lag_flag[i] = 0 
                lag_begins[i] = 0 # there might be multiple bugs doing their lags together, so need to screen everything again
    if(ignore_this_run == False):
        global details
        details['t_info'].append(t_points)
        details['bug_info'].append((survivors, c_points))
        details['res_left'].append(r_points)
        if('res_concent' in details):
            details['res_concent'].append(conc_points)
        if('round_idx' in details.keys()):
            details['round_idx'].append(details['round_idx'][-1]+1)
    return system, t_points, c_points, r_points
# you can modify the fluctuating environment here.
def move_to_new(system, fluc = []):
    global Res
    global details
    global movestep
    movestep += 1
    baseline = [Nr * i / sum(Res) for i in Res]
    system['bug_density'] = [(i*D > b_threshold)*i*D for i in system['bug_density']]
    system['bug_available'] = [1*(system['bug_density'][i] > b_threshold) for i in range(Nb)]
    system['res_concentration'] = [(system['res_concentration'][i]*D + baseline[i]) / (D + 1) for i in range(Nr)]
    system['res_available'] = [1*(system['res_concentration'][i] > c_threshold) for i in range(Nr)]
    Res = baseline
    details['res_begin'].append(baseline)
    return system

def invade(system, bug, growth_rate_list, invlist):
    global details
    # (starting from a new flask)
    # then add invasive species
    system1 = copy.deepcopy(system)
    # then dilute till steady
    system1, t_points, c_points, r_points = dilute_check(system1, growth_rate_list, lagTime, True)
    #system1, t_points, c_points, r_points = dilute_check(system1, growth_rate_list)
    accum = 0
    if(len(t_points)>2):
        t_begin = t_points[0]
        i = 0
        nut = 0
        for idx in range(len(t_points[:-1])):
            switch_flag = 0 # if switch to next nutrient
            t_end = t_points[idx + 1]
            pref_bug = preference_list[bug//Size]
            while(r_points[idx][pref_bug[i]] <= 0 and i < Nr-1): 
                i += 1
                switch_flag = 1
            if(idx == len(t_points[:-1]) - 1):
                switch_flag = 1
            if(switch_flag == 1):
                if(t_begin == t_points[0]):
                    t_slot = t_end - t_begin
                else:
                    t_slot = t_end - t_begin - lagTime
                if(t_slot > 0):
                    accum += t_slot*growth_rate_list[bug][preference_list[bug//Size][nut]]
                t_begin = t_end
                nut += 1
    else:
        t_slot = t_points[1]-t_points[0]
        accum = t_slot*growth_rate_list[bug][preference_list[bug//Size][0]]
    if accum > log(1/D) or sum(system['bug_available']) == 0:
        system['bug_available'][bug] = 1
        system['bug_density'][bug] = b0
        i=0
        eqm = 0
        while(len(output(system)[0])>Nr or i < dilute_to_steady or eqm == 0):
            i=i+1
            system, t_points, c_points, r_points = dilute_check(system, growth_rate_list)
            if sum(system['bug_available']) == 0:
                details = {'res_begin':[], 't_info':[], 'res_left':[], 'bug_info':[], 'round_idx':[0]}
                system = {'res_available': np.heaviside(Res, 0), 'res_concentration': [i for i in Res], 'bug_available': [0 for i in range(Nb)], 'bug_density': [0 for i in range(Nb)]}
                break
            # move to a new flask
            system = move_to_new(system, fluc=[])
            if(len(details['bug_info']) > 2):
                if(len(details['bug_info'][-1][0]) == len(details['bug_info'][-2][0])):
                    if max([abs(details['bug_info'][-1][1][1][i]-details['bug_info'][-2][1][1][i])/details['bug_info'][-2][1][1][i] for i in range(len(details['bug_info'][-2][1][1]))]) < 5e-3:
                        eqm = 1
        system, t_points, c_points, r_points = dilute_check(system, growth_rate_list)
        system = move_to_new(system, fluc=[])
        ext_list = [i for i in invlist if system['bug_available'][i] == 0]
        return system, ext_list
    else:
        ext_list = [i for i in invlist if system['bug_available'][i] == 0]
        if('round_idx' in details.keys()):
            details['round_idx'][-1] += T_dilute
        return system, ext_list

def round_robin_invade(system, ext_list, growth_rate_list, lagTime):
    for bug in ext_list:
        system, new_ext_list = invade(system, bug, growth_rate_list, ext_list)
    return system, new_ext_list

def output(system):
    survivors = [i for i, v in enumerate(system['bug_available']) if v != 0]
    pref_list = []
    concent = []
    growth = []
    for i in survivors:
        pref_list.append(preference_list[i//Size])
        concent.append(system['bug_density'][i])
        growth.append(growth_rate_list[i])
    return survivors, concent, pref_list, growth

def MeasureOutput(output, resources, g):
    system = {'res_available': np.heaviside(resources, 0), 'res_concentration': [i for i in resources], 'bug_available': [0 for i in range(Nb)], 'bug_density': [0 for i in range(Nb)]}
    for bug in output[0]:
        system['bug_available'][bug] = 1
        system['bug_density'][bug] = output[1][output[0].index(bug)]
    for i in range(1):
        system, p1, p2, p3 = dilute_check(system, g)
        baseline = resources
        baseline = [Nr * i / sum(baseline) for i in baseline]
        system['bug_density'] = [(i*D > b_threshold)*i*D for i in system['bug_density']]
        system['bug_available'] = [1*(system['bug_density'][i] > b_threshold) for i in range(Nb)]
        system['res_concentration'] = [(system['res_concentration'][i]*D + baseline[i]) / (D + 1) for i in range(Nr)]
        system['res_available'] = [1*(system['res_concentration'][i] > c_threshold) for i in range(Nr)]
    system, p1, p2, p3 = dilute_check(system, g)
    return p1, p2, p3

def get_experts(growth_rate_list):
    Nb = len(growth_rate_list)
    Nr = len(growth_rate_list[0])
    experts = [[0, 0] for i in range(Nr)]
    for bug in range(len(growth_rate_list)):
        if(growth_rate_list[bug][int(bug//(Nb/Nr))]>experts[int(bug//(Nb/Nr))][1]):
            experts[int(bug//(Nb/Nr))][1] = growth_rate_list[bug][int(bug//(Nb/Nr))]
            experts[int(bug//(Nb/Nr))][0] = bug
    return experts

def crit1(e):
    return e[1]
def crit0(e):
    return e[0]
def get_alternate_experts(growth_rate_list):
    Nb = len(growth_rate_list)
    Nr = len(growth_rate_list[0])
    topChoiceList = []
    for bug in range(len(growth_rate_list)):
        topChoiceList.append((bug, growth_rate_list[bug][int(bug//(Nb/Nr))]))
    topChoiceList.sort(key=crit1, reverse = True)
    experts = topChoiceList[:Nr]
    experts.sort(key=crit0)
    return experts

def isFreak(growth_rate_list, bug):
    Nb = len(growth_rate_list)
    Nr = len(growth_rate_list[0])
    if growth_rate_list[bug][int(bug//(Nb/Nr))] != max(growth_rate_list[bug]):
        return True
    return False

def FreaksFromOutput(output):
    freaks = []
    for i in range(len(output[0])):
        PrefOrder = output[2][i]
        g = output[3][i]
        if(max(g)!=g[PrefOrder[0]]):
            freaks.append(output[0][i])
    return freaks

def FreakGrowthRates(output):
    freakGs = []
    for i in range(len(output[0])):
        PrefOrder = output[2][i]
        g = output[3][i]
        if(max(g)!=g[PrefOrder[0]]):
            freakGs.append(output[3][i])
    return freakGs

# generating growth rates of different distributions

def G(Nb, Nr, budget, g_var):
    growth_rate_list = np.random.random([Nb, Nr])*np.exp(np.random.normal(0, g_var, (Nb, 1)))*budget
    filt = np.heaviside(growth_rate_list, 0)
    growth_rate_list = filt*growth_rate_list
    #growth_rate_list = budget*growth_rate_list / np.sum(growth_rate_list, axis=1)[:, None]*np.exp(np.random.normal(0, g_var, (Nb, 1)))
    return growth_rate_list

def GLN(Nb, Nr, budget, g_var):
    growth_rate_list = np.exp(np.random.normal(0, g_var, (Nb, Nr)))*budget
    return growth_rate_list


def GFreak(Nb, Nr, budget, g_var, bias):
    growth_rate_list = np.random.random([Nb, Nr])
    filt = np.heaviside(growth_rate_list, 0)
    growth_rate_list = filt*growth_rate_list
    growth_rate_list = budget*growth_rate_list / np.sum(growth_rate_list, axis=1)[:, None]*np.exp(np.random.normal(0, g_var, (Nb, 1)))
    for i in range(len(growth_rate_list)):
        first = int(i//(Size*math.factorial(Nr)/Nr))
        if growth_rate_list[i][first] != max(growth_rate_list[i]):
            growth_rate_list[i] = [elem + bias for elem in growth_rate_list[i]]
    return growth_rate_list

def MakeGVeryNonFreak(G):
    GOut = []
    for g in G:
        gout = list(g)
        gout.sort(reverse = True)
        GOut.append(gout)
    return GOut

def PickInvadersFor2Res(G):
    InvList = []
    NPrefOrder = 2
    Nb = len(G)
    Pref1 = G[:int((Nb+1)/2)]
    expert1 = get_experts(G)[0][0]
    InvList.append(expert1)
    for idx, i in enumerate(Pref1):
        if i[1] > G[expert1][1]:
            InvList.append(idx)
    Pref2 = G[int((Nb+1)/2):]
    expert2 = get_experts(G)[1][0]
    InvList.append(expert2)
    for idx, i in enumerate(Pref2):
        if i[0] > G[expert2][0]:
            InvList.append(idx+int((Nb+1)/2))
    return InvList

# about anomalous species
def tellFreaks2d(species_list, Size, v = False):
    freaks_list = []
    for species in species_list:
        if species < Size:
            #then it consumes nut 0 first
            if growth_rate_list[species][1]>growth_rate_list[species][0]:
                freaks_list.append(species)
        else:
            #then it consumes nut 1 first
            if growth_rate_list[species][1]<growth_rate_list[species][0]: freaks_list.append(species)
    if v == True:
        alter_growth_rate_list = growth_rate_list.copy()
        for species in species_list:
            if species >= Size: alter_growth_rate_list[species][0], alter_growth_rate_list[species][1] = alter_growth_rate_list[species][1], alter_growth_rate_list[species][0]
        plt.scatter([alter_growth_rate_list[species][0] for species in species_list], [alter_growth_rate_list[species][1] for species in species_list], color = 'b')
        plt.scatter([alter_growth_rate_list[species][0] for species in freaks_list], [alter_growth_rate_list[species][1] for species in freaks_list], color = 'b')
        plt.plot(np.arange(0, 9, 0.01), np.arange(0, 9, 0.01))
        plt.xlabel('top choice g')
        plt.ylabel('second choice g')
    return freaks_list

def makeRegular(G):
    Nb = len(G)
    for g in G[:int(Nb/2)]:
        g[0], g[1] = max(g[0], g[1]), min(g[0], g[1])
    for g in G[int(Nb/2):]:
        g[1], g[0] = max(g[0], g[1]), min(g[0], g[1])
    return G

# about complementarity
def preproc(outputs):
    S = []
    for i in range(len(outputs[0][2][0])):
        S_data = []
        for output_case in outputs:
            S_data.append([preflist[i] for preflist in output_case[2]])
        S.append(S_data)
    return S

def getEntropy(list0):
    return len(list(set(list0)))/len(list0)

def getEntropies(S, R):
    avr = []
    std = []
    for elem in S:
        avr.append(np.mean([getEntropy(lst) for lst in elem if len(lst) == R]))
        std.append(sqrt(np.var([getEntropy(lst) for lst in elem if len(lst) == R]))/sqrt(len(elem)))
    rdm = np.mean([(1-(1-1/R)**len(lst))*R/len(lst) for lst in S[0] if len(lst) == R])
    return avr, std, rdm

# efficient market
def TPoints2Taus(t_points):
    Taus = []
    for idx in range(len(t_points)-2):
        Taus.append(t_points[idx+1]-t_points[idx])
    return Taus

def RPoints2DepOrder(r_points):
    DepOrder = []
    for rpoint in range(len(r_points[1:])):
        subtracted = [int(r_points[rpoint][i] - r_points[rpoint+1][i]) for i in range(len(r_points[0]))]
        if 1 not in subtracted:
            return []
        DepOrder.append(subtracted.index(1))
    return DepOrder

# get the necessary informations along the assembly process
def getMarketHistory(details):
    # first list is a marker
    # last list is the real time where this invasion is supposed to happen. In our invasion codes, uninvasible bugs are never allowed to get in, but in reality they take some time to fail the invasions.
    # print(details['bug_info'])
    bugstates = [[1], [details['bug_info'][0][0]], [RPoints2DepOrder(details['res_left'][0])], [TPoints2Taus(details['t_info'][0])], [details['round_idx'][0]], [details['bug_info'][0][1]]]
    for idx, elements in enumerate(details['bug_info'][1:]):
        if(elements[0]!=bugstates[1][-1]):
            bugstates[0].append(1)
            bugstates[1].append(elements[0])
            bugstates[2].append(RPoints2DepOrder(details['res_left'][idx]))
            bugstates[3].append(TPoints2Taus(details['t_info'][idx]))
            bugstates[4].append(details['round_idx'][idx])
            bugstates[5].append(elements[1])
        else:
            bugstates[0][-1]+=1
            bugstates[1][-1] = (elements[0])
            bugstates[2][-1] = (RPoints2DepOrder(details['res_left'][idx]))
            bugstates[3][-1] = (TPoints2Taus(details['t_info'][idx]))
            bugstates[4][-1] = details['round_idx'][idx]
            bugstates[5][-1] = (elements[1])
    # and make sure all are equilibrium. first, eliminate all the mayflies (gets in for only a short time,
    # because using g*t vs logD to tell "if gets in" has errors: we assumed B0=0, where actually it's b0)
    for idx, i in enumerate(bugstates[0]):
        if i < 2:
            del bugstates[0][idx]
            del bugstates[1][idx]
            del bugstates[2][idx]
            del bugstates[3][idx]
            del bugstates[4][idx]
            del bugstates[5][idx]
    # and second, eliminate the intermediate stages, where bugs take time to die out but no one is being introduced into system
    idx = 0
    while(idx < len(bugstates[1])-1):
        if set(bugstates[1][idx])>=set(bugstates[1][idx+1]):
            del bugstates[0][idx]
            del bugstates[1][idx]
            del bugstates[2][idx]
            del bugstates[3][idx]
            del bugstates[4][idx]
            del bugstates[5][idx]
            idx-=1
        idx+=1
    idx = 0
    while(idx < len(bugstates[1])-1):
        if set(bugstates[1][idx])>=set(bugstates[1][idx+1]):
            del bugstates[0][idx]
            del bugstates[1][idx]
            del bugstates[2][idx]
            del bugstates[3][idx]
            del bugstates[4][idx]
            del bugstates[5][idx]
            idx-=1
        idx+=1
    # and finally, list the invaders
    bugstates.append([bugstates[1][0][0]])
    for idx, i in enumerate(bugstates[1][1:]):
        for j in i:
            if j not in bugstates[1][idx-1+1]: # idx still starts from 0
                bugstates[-1].append(j)
    return bugstates

def plotTimes(data):
    ax = plt.gca()
    for k in range(len(data[2][0])):
        ax.plot(range(len(data[3])), [i[k] for i in data[3]])
    ax.plot(range(len(data[3])), [sum(i) for i in data[3]], label = "overall time")
    ax.set_xticklabels(range(len(data[2])+1))
    plt.xlabel("invasion")
    plt.ylabel("tau")
    plt.legend()

def PrefList2ConsumeList(plist, dlist):
    preflist = [k for k in plist]
    depletelist = [k for k in dlist]
    consumelist = []
    while(preflist!=[]):
        if(preflist[0] in depletelist):
            consumelist.append(preflist[0])
            del depletelist[0]
        else:
            del preflist[0]
    return consumelist