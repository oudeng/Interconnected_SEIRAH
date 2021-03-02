#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb., 21, 2021
@author:
    (Structure & Algorithms) Deng
    (Coding) Deng, Sandy

    MIT License
    Copyright (c) 2021 Ou, DENG (dengou@toki.waseda.jp)

Name:
    Main of Interconnected SEIRAH Framework

Purpose:
    Experiment attached to research paper titled “An Extended Epidemic Model
    on Interconnected Networks for COVID-19 to Explore the Epidemic Dynamics”. 
    The paper was coauthored by Ou Deng, Kiichi Tago, and Qun Jin. 
    (paper1 ver1.9.2)

Environment:
    1. Software: Python 3.7 with NetworkX library.
    2. Hardware: Alicloud (high perfomance computing is necessary)
    
Note:
    Rt, Tt need the larger network than current experiment.
    
"""

import networkx as nx
import random
from random import sample
import matplotlib.pyplot as plt1
import matplotlib.pyplot as plt2
import matplotlib.pyplot as plt3
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import copy as co
import csv as cs
import pickle as pic
import datetime
import SEIRAH_SW as SEIRAH_SW
import SEIRAH_SW_F as SEIRAH_SW_F

random.seed(2020)

"""
Network Generation and Initialization.
The newman_watts_strogatz_graph as following official document.
https://networkx.github.io/documentation/stable/_modules/networkx/generators/random_graphs.html#newman_watts_strogatz_graph
"""
def G_gene(N,k,p,expo,infe,asym,hosp,reco):
    G = nx.newman_watts_strogatz_graph(N,k,p,seed=2020);  
    
    for i in G.nodes():
        G.nodes[i]['status'] = 'susc'  #Initialize all nodes as susc.
        G.nodes[i]['S_1stday'] = 0  #Set initial date(day) to S nodes. 
        G.nodes[i]['E_1stday'] = 0
        G.nodes[i]['I_1stday'] = 0
        G.nodes[i]['A_1stday'] = 0
        G.nodes[i]['H_1stday'] = 0
        G.nodes[i]['R_1stday'] = 0
        G.nodes[i]['Infe_other'] = 0  #Count for Rt calculating. 
        
    for i in G.edges():
        G.edges[i]['weight'] = 1   #Initialize all edges weight = 1 
    
    """
    Initial G. SEIRAH status
    """
    samp_sum = expo+infe+asym+hosp+reco   #total number of samples

    
    # Random sample nodes from G.
    random_n = sample(G.nodes(), samp_sum)  #use sample function
    
    # nodes lists of each status. Randomly chose(clip) from random_n list.
    node_e = random_n[0:expo]
    
    a = expo + infe
    node_i = random_n[expo:a]
    
    b = a + asym
    node_a = random_n[a:b]
    
    c = b + hosp
    node_h = random_n[b:c]
    
    d = c + reco
    node_r = random_n[c:d]
    
    # Divide samples by status inputed numbers. Assign status except S.
    # Record *_1stday. S node does not need to record.
    for i in node_e:
        G.nodes[i]['status'] = 'expo'
        G.nodes[i]['S_1stday'] = 0 #Initial E to replace S status.
        G.nodes[i]['E_1stday'] = 0
        G.nodes[i]['I_1stday'] = 0
        G.nodes[i]['A_1stday'] = 0
        G.nodes[i]['H_1stday'] = 0
        G.nodes[i]['R_1stday'] = 0
    
    for i in node_i:
        G.nodes[i]['status'] = 'infe'
        G.nodes[i]['S_1stday'] = 0 #Initial E to replace S status.
        G.nodes[i]['E_1stday'] = 0
        G.nodes[i]['I_1stday'] = 0
        G.nodes[i]['A_1stday'] = 0
        G.nodes[i]['H_1stday'] = 0
        G.nodes[i]['R_1stday'] = 0
        
    for i in node_a:
        G.nodes[i]['status'] = 'asym'
        G.nodes[i]['S_1stday'] = 0 #Initial E to replace S status.
        G.nodes[i]['E_1stday'] = 0
        G.nodes[i]['I_1stday'] = 0
        G.nodes[i]['A_1stday'] = 0
        G.nodes[i]['H_1stday'] = 0
        G.nodes[i]['R_1stday'] = 0
        
    for i in node_h:
        G.nodes[i]['status'] = 'hosp'
        G.nodes[i]['S_1stday'] = 0 #Initial E to replace S status.
        G.nodes[i]['E_1stday'] = 0
        G.nodes[i]['I_1stday'] = 0
        G.nodes[i]['A_1stday'] = 0
        G.nodes[i]['H_1stday'] = 0
        G.nodes[i]['R_1stday'] = 0
    
    for i in node_r:
        G.nodes[i]['status'] = 'reco'
        G.nodes[i]['S_1stday'] = 0 #Initial E to replace S status.
        G.nodes[i]['E_1stday'] = 0
        G.nodes[i]['I_1stday'] = 0
        G.nodes[i]['A_1stday'] = 0
        G.nodes[i]['H_1stday'] = 0
        G.nodes[i]['R_1stday'] = 0
        
    # update G edges weight. node_h needs update link weight 1=>0
    for i in node_h:
        for nbr in G[i]:
            G.edges[i, nbr]['weight'] = 0            

    return G  

"""
Nodes mapping function. Shift weight as 0 or 1.
"""
def edge_weight_0(nwk):
    for i in nwk.nodes():       
        if nwk.nodes[i]['status'] == 'hosp' or 'Com' in nwk.nodes[i]:
            for nbr in nwk[i]:
                nwk.edges[i, nbr]['weight'] = 0          
    return nwk

def edge_weight_1(nwk):
    for i in nwk.nodes():       
        if nwk.nodes[i]['status'] != 'hosp' or 'Com' in nwk.nodes[i]:           
            for nbr in nwk[i]:
                nwk.edges[i, nbr]['weight'] = 1          
    return nwk

"""
Search for optimal beta_t. Step of Algorithm2.
"""
def find_beta_t(__day,__beta_t,__CBD_w,__city_0,__city_1,__city_2,__city_3,
                __cityComList_0,__cityComList_1,__cityComList_2,__cityComList_3):
    
    eps=0.01  #episilon in Algorithm2
    beta_p=__beta_t
    beta_p_b=0
    
    df = pd.read_csv('new_cases_cr2020.csv', nrows=day_n+8)
    days = df[__day:(__day+7)]
    
    real_H = np.zeros((7,8))
    
    real_H[:,5]=days['H_Shuto']
    
    __CBD_w_=pic.loads(pic.dumps(__CBD_w))
    __city_0_=pic.loads(pic.dumps(__city_0))
    __city_1_=pic.loads(pic.dumps(__city_1))
    __city_2_=pic.loads(pic.dumps(__city_2))
    __city_3_=pic.loads(pic.dumps(__city_3))
    __cityComList_0_=pic.loads(pic.dumps(__cityComList_0))
    __cityComList_1_=pic.loads(pic.dumps(__cityComList_1))
    __cityComList_2_=pic.loads(pic.dumps(__cityComList_2))
    __cityComList_3_=pic.loads(pic.dumps(__cityComList_3))
    
    A = 0  #beta_t in [0,1]
    B = 1

    try:
        while abs(beta_p - beta_p_b) >= eps : # 
        
            #print(A,beta_p,B)
        
            beta_p = (A+B)/2
            
            beta_t_l=[A,beta_p,B]
            
            #print(A,beta_p,B,beta_p_b)
            
            __CBD_w=pic.loads(pic.dumps(__CBD_w_))
            __city_0=pic.loads(pic.dumps(__city_0_))
            __city_1=pic.loads(pic.dumps(__city_1_))
            __city_2=pic.loads(pic.dumps(__city_2_))
            __city_3=pic.loads(pic.dumps(__city_3_))
            __cityComList_0=pic.loads(pic.dumps(__cityComList_0_))
            __cityComList_1=pic.loads(pic.dumps(__cityComList_1_))
            __cityComList_2=pic.loads(pic.dumps(__cityComList_2_))
            __cityComList_3=pic.loads(pic.dumps(__cityComList_3_))
               
            D_A=beta_leastsquare(beta_t_l[0],__CBD_w,__city_0,__city_1,__city_2,__city_3,
                                    __cityComList_0,__cityComList_1,__cityComList_2,__cityComList_3,
                                    real_H)
            
            __CBD_w=pic.loads(pic.dumps(__CBD_w_))
            __city_0=pic.loads(pic.dumps(__city_0_))
            __city_1=pic.loads(pic.dumps(__city_1_))
            __city_2=pic.loads(pic.dumps(__city_2_))
            __city_3=pic.loads(pic.dumps(__city_3_))
            __cityComList_0=pic.loads(pic.dumps(__cityComList_0_))
            __cityComList_1=pic.loads(pic.dumps(__cityComList_1_))
            __cityComList_2=pic.loads(pic.dumps(__cityComList_2_))
            __cityComList_3=pic.loads(pic.dumps(__cityComList_3_))
            
            D_B=beta_leastsquare(beta_t_l[2],__CBD_w,__city_0,__city_1,__city_2,__city_3,
                                    __cityComList_0,__cityComList_1,__cityComList_2,__cityComList_3,
                                    real_H)
            
            #print(D_A,D_B)
            
            if D_A > D_B:
                A = co.copy(beta_t_l[1])
            if D_A < D_B:
                B = co.copy(beta_t_l[1])
            
            beta_p_b= (A+B)/2
            
            #print(A,B)
            
            #print(beta_p,beta_p_b)
    
    except Exception as re:
        print("Terminated by abnormal network topology")
        print("Error Message:",re)
    
    return beta_p

"""
Simulate β by H(t1,t2,...t7). MSE method.
"""

def beta_leastsquare(beta_t,__CBD_s,__city_0,__city_1,__city_2,__city_3,
                     __cityComList_0,__cityComList_1,__cityComList_2,__cityComList_3,
                     __real_H):
    
    __Daily_Result_Week=np.array([0,0,0,0,0,0,0,0])
    
    CBD_s=pic.loads(pic.dumps(__CBD_s))    
    
    vars()["citys_" + str(0)]=pic.loads(pic.dumps(__city_0))
    vars()["citys_" + str(1)]=pic.loads(pic.dumps(__city_1))
    vars()["citys_" + str(2)]=pic.loads(pic.dumps(__city_2))
    vars()["citys_" + str(3)]=pic.loads(pic.dumps(__city_3))
    
    vars()["cityComLists_" + str(0)]=pic.loads(pic.dumps(__cityComList_0))
    vars()["cityComLists_" + str(1)]=pic.loads(pic.dumps(__cityComList_1))
    vars()["cityComLists_" + str(2)]=pic.loads(pic.dumps(__cityComList_2))
    vars()["cityComLists_" + str(3)]=pic.loads(pic.dumps(__cityComList_3))

    for d in range(7):  # Day: from 0 to 6
        for i in range(4):  #all 4 cities in tau[0] TimeZone
            __city_for_func = vars()["citys_" + str(i)] # prepare which city to do.
            vars()["citys_" + str(i)] = SEIRAH_SW_F.SEIRAH_SW_F(__city_for_func,beta_t,tau[0])
            
            # Finish tau[0]:life infection. Prepare for tau[1]:life+CBD
            edge_weight_0(__city_for_func)
            
        #Interconnect_City2CBD()  
        count = 0
        for i in range(4):
            for j in vars()["cityComLists_" + str(i)]:
                __status_value = vars()["citys_" + str(i)].nodes[j]['status']
                CBD_s.nodes[count]['status'] = __status_value
                count = count + 1 
        
        # SEIRAH process of CBD in a certain day.
        CBD_s = SEIRAH_SW_F.SEIRAH_SW_F(CBD_s,beta_t, tau[1])
        
        # SEIRAH precess of cities in working TimeZone of a certain day.    
        for i in range(4):  #all 4 cities in TimeZone tau[1]
            __city_for_func = vars()["citys_" + str(i)] # prepare which city to do.
            vars()["citys_" + str(i)] = SEIRAH_SW_F.SEIRAH_SW_F(__city_for_func,beta_t,tau[1])
            edge_weight_1(__city_for_func)
            
        #Interconnect_CBD2City()      
        count = 0
        for i in range(4):
            for j in vars()["cityComLists_" + str(i)]:
                __status_value = CBD_s.nodes[count]['status']
                vars()["citys_" + str(i)].nodes[j]['status'] = __status_value
                count = count + 1
                
        #Output daily status. Stoped S counting.
        __Daily_Result_sum = np.array([0,0,0,0,0,0,0,0]) # S, E, I, R, A, H, Rt, Tt
        
        for i in range(4):
            __city_for_func = vars()["citys_" + str(i)] 
            __Daily_Result_city = Count_status(day,__city_for_func)
            __Daily_Result_sum = __Daily_Result_sum + __Daily_Result_city    
            
        
        if d == 0:
            __Daily_Result_Week =__Daily_Result_sum
        else:
            __Daily_Result_Week = np.vstack((__Daily_Result_Week, __Daily_Result_sum))
        
    
    __Daily_Result_df = pd.DataFrame(__Daily_Result_Week)
        
    pred_H = np.zeros((7,8))
    pred_H[:,5] = __Daily_Result_df.iloc[:,5]
    
    diff=(__real_H-pred_H)**2
        
    sum_diff=0
    
    for i in range(1,7): 
        sum_diff = sum_diff + diff[i,5]
    
    sum_diff_s=sum_diff
    
    return sum_diff_s

def Count_status(__time_stamp, __nwk):
    
    num_S = num_E = num_I = num_R = num_A = num_H = 0
    count_Infe_other = count_T = 0
    Rt = Tt = 0

    for i in range(__nwk.number_of_nodes()):
# S is big number. No count to save computing time.
#        if nwk.nodes[i]['status'] == 'susc':
#            num_S = num_S + 1
#           continue
        
        if __nwk.nodes[i]['status'] == 'expo':
            num_E = num_E + 1
            #continue
        if __nwk.nodes[i]['status'] == 'infe':
            num_I = num_I + 1
            #continue
        if __nwk.nodes[i]['status'] == 'reco':
            num_R = num_R + 1
            #continue
        if __nwk.nodes[i]['status'] == 'asym':
            num_A = num_A + 1
            #continue
        if __nwk.nodes[i]['status'] == 'hosp':
            num_H = num_H + 1
            
        #Filter out new H, to count Rt, Tt
        if __nwk.nodes[i]['H_1stday'] == __time_stamp:
            count_Infe_other = count_Infe_other + __nwk.nodes[i]['Infe_other']
            #Count how many days from E to H
            count_T = count_T + __nwk.nodes[i]['H_1stday'] - __nwk.nodes[i]['E_1stday']

        if num_H == 0:  #Avoid error while H=0
            Rt = Tt = 0
    
        else:
            Rt = count_Infe_other / num_H
            Tt = count_T / num_H
            #Tt = format(Rt / beta_t,'.2f')
    
    __Count_status = [num_S, num_E, num_I, num_R, num_A, num_H, Rt, Tt]
    
    return __Count_status


# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #       Above sections are FUNCTIONS.             # #
# #       MAIN() is as follows.                     # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

'''
Spatiotemporal Analysis of SEIRAH on interconnected small world.
'''
# Parameter input section. Global variables for functions use easily.
# More important, SEIRAH_SW and SEIRAH use the same parameters settings.
# Parameter meaning:
# N:population. k:SW node degree. p:rewired probability. <=SW parameters.
# NCom: number of commuting people.

para_0 = {'N':135200,'k':4,'p':0.005,'NCom':48640,
          'expo':3,'infe':0,'asym':0,'hosp':0,'reco':0}
    
para_1 = {'N':92000,'k':4,'p':0.005,'NCom':8880,
          'expo':3,'infe':0,'asym':0,'hosp':0,'reco':0}
    
para_2 = {'N':73400,'k':4,'p':0.005,'NCom':7800,
          'expo':3,'infe':0,'asym':0,'hosp':0,'reco':0}

para_3 = {'N':62800,'k':4,'p':0.005,'NCom':5980,
          'expo':3,'infe':0,'asym':0,'hosp':0,'reco':0}

"""
# Use following small N to test. Big N cost too much computing time.
para_0 = {'N':1352,'k':4,'p':0.1,'NCom':486,
          'expo':3,'infe':0,'asym':0,'hosp':0,'reco':0}
    
para_1 = {'N':920,'k':4,'p':0.1,'NCom':88,
          'expo':3,'infe':0,'asym':0,'hosp':0,'reco':0}
    
para_2 = {'N':734,'k':4,'p':0.1,'NCom':78,
          'expo':3,'infe':0,'asym':0,'hosp':0,'reco':0}

para_3 = {'N':628,'k':4,'p':0.1,'NCom':59,
          'expo':3,'infe':0,'asym':0,'hosp':0,'reco':0}
"""

para_N_0=co.copy(para_0)
para_N_1=co.copy(para_1)
para_N_2=co.copy(para_2)
para_N_3=co.copy(para_3)

para_Epidemic = {'sigma':0.2,'p1':0.18,'p2':0.3,'l_AH':0.05,'l_IH':0.3,
                 'g_AR':0.07,'g_HR':0.1}

#beta = [0.1,0.1,0.1] # B_IS, B_AS, B_ES for future upgrade

tau = [0.5, 0.5] # Default 0.2,0.4 TimeZone(1), TimeZone(2)...Rest TimeZone Not actively infection

# (Epidemic parameters) ###################################################
# beta: Latency rate. β* in equations.
# sigma: Transmission rate. σ in equations.
# p1: Ratio of E to A.
# l_*: Hospitalization rate. λ* in equations.
# p2: Hospitalization ratio of A.
# g_*: Recovery rate. γ* in equations.
# tau: TimeZone Ratio Coefficient. τ* in equations.

# CBD is Commuting network. Interconnnected with city_0,1,2,3
CBD_N = para_0['NCom']+para_1['NCom']+para_2['NCom']+para_3['NCom']
CBD_k = 8   # default k=6
CBD_p = 0.05  # default CBD_p = 0.5

cityComList =[]
l=0

# Initial day, β, by Sandy
day_n=60 #Initial days
arr_beta = np.arange(day_n)
beta_t_S = pd.Series(0.1,arr_beta)
## parameters input END ######

'''
Initialize commuting network
'''
#Commuting network. nx.draw(H,pos=nx.circular_layout(H))
CBD = nx.newman_watts_strogatz_graph(CBD_N,CBD_k,CBD_p,seed=2980) 
for i in CBD.nodes():
    CBD.nodes[i]['Infe_other'] = 0  #Count for Rt calculating. 

for i in CBD.edges():
    CBD.edges[i]['weight'] = 1


'''
Use vars() to transfer string to variable.
Network of city_0,1,2....: city_0 is center, the others are outskirts. 
'''
count = 0  # CBD nodes index.
for i in range(4):  # 4 cities.
    
    # Use parameter set of target city.
    para = vars()["para_" + str(i)]
    para_N = vars()["para_N_" + str(i)]
    
    #print(para_N)
    
    vars()["city_" + str(i)] = G_gene(para['N'],para['k'],para['p'],para['expo'],
                                      para['infe'],para['asym'],para['hosp'],para['reco'])
    '''
    Choose cityCom（通勤者）nodes. Output lists of cityComList_0,1,2,3
    '''
    # From city_* to choose sample nodes, with sample number='NCom'.  
    # Later Ncom will be daily changed by CR(commuting ratio)
    vars()["cityComList_" + str(i)] = sample(vars()["city_" + str(i)].nodes(),para['NCom'])

    # cityCom is all lists of commuting nodes of cityCom_*
    cityComList = cityComList + vars()["cityComList_" + str(i)]
    N_cityComList = len(cityComList)  # = NH
    
    # cityCom_* is nodes ID list of commuting. Add 'Com':True to those nodes in city_* 
    for j in vars()["cityComList_" + str(i)]:
        vars()["city_" + str(i)].nodes[j]['Com'] = True
        # Copy commting nodes 'status' to CBD nodes. 
        status_value = vars()["city_" + str(i)].nodes[j]['status']
        CBD.nodes[count]['status'] = status_value
        count = count + 1

# While generated, all CBD.edges as 'weight'=1
# Initialize 'hosp' status nodes with nbr edges of 'weight'=0        
edge_weight_0(CBD)

'''
Above finished Networks initialization for city_* and CBD.
SEIRAH_SW and SEIRAH will use the same initial parameter settings.
'''

'''
Now Start Spatiotempral SEIRAH_SW Status Transition.
'''
#Initial Status
day = 0
Daily_Result_sum = np.array([0,0,0,0,0,0,0,0]) 
# Initial S,E,I,R,A,H,Rt,Tt

tmp_S = tmp_E = tmp_I = tmp_R = tmp_A = tmp_H = 0

for i in range(4):
    city_for_func = vars()["city_" + str(i)] # prepare which city to do.
    Daily_Result_city = Count_status(day,city_for_func)
    print('city_', i, Daily_Result_city)
    Daily_Result_sum = Daily_Result_sum + Daily_Result_city

print('day:', day, 'Total',Daily_Result_sum)

Daily_Result_ALL = Daily_Result_sum

cityresult_0= np.array([0,0,0,0,0,0,0,0]) 
cityresult_1= np.array([0,0,0,0,0,0,0,0])
cityresult_2= np.array([0,0,0,0,0,0,0,0])
cityresult_3= np.array([0,0,0,0,0,0,0,0])

print()

df = pd.read_csv('new_cases_cr2020.csv', nrows=day_n+8)

Ncom_l = df['CR_Shuto']

print('Start:')

for day in range(day_n):  # Day: from 0 to (n-1)

    count = 0  # CBD nodes index.
    
    for cr in range(4):  # 4 cities.
    
        para = vars()["para_" + str(cr)]
        para_N = vars()["para_N_" + str(cr)]
        
        #print(para,para_N)

        para['NCom'] = int(para_N['NCom']*Ncom_l[day])
        
        #print(para['N'],Ncom_l[day],para['NCom'])

        vars()["cityComList_" + str(cr)] = sample(vars()["city_" + str(cr)].nodes(),para['NCom'])
        
        cityComList = cityComList + vars()["cityComList_" + str(cr)]
        N_cityComList = len(cityComList)  # = NH

        # cityCom_* is nodes ID list of commuting. Add 'Com':True to those nodes in city_* 
        for j in vars()["cityComList_" + str(cr)]:
            vars()["city_" + str(cr)].nodes[j]['Com'] = True
            # Copy commting nodes 'status' to CBD nodes. 
            status_value = vars()["city_" + str(cr)].nodes[j]['status']
            #print(status_value)
            CBD.nodes[count]['status'] = status_value
            count = count + 1
    
    edge_weight_0(CBD)

    if day ==0:
        beta_t=beta_t_S.iloc[day]
        for i in range(4):  #all 4 cities in tau[0] TimeZone
            city_for_func = vars()["city_" + str(i)] # prepare which city to do.
            SEIRAH_SW.SEIRAH_SW(day,city_for_func,beta_t,tau[0])
            
            # Finish tau[0]:life infection. Prepare for tau[1]:life+CBD
            edge_weight_0(city_for_func)
            
        #Interconnect_City2CBD()  
        count = 0
        for i in range(4):
            for j in vars()["cityComList_" + str(i)]:
                status_value = vars()["city_" + str(i)].nodes[j]['status']
                CBD.nodes[count]['status'] = status_value
                count = count + 1 
        
        # SEIRAH process of CBD in a certain day.
        SEIRAH_SW.SEIRAH_SW(day,CBD,beta_t, tau[1])
        
        # SEIRAH precess of cities in working TimeZone of a certain day.    
        for i in range(4):  #all 4 cities in TimeZone tau[1]
            city_for_func = vars()["city_" + str(i)] # prepare which city to do.
            SEIRAH_SW.SEIRAH_SW(day,city_for_func,beta_t,tau[1])
            edge_weight_1(city_for_func)
            
        #Interconnect_CBD2City()      
        count = 0
        for i in range(4):
            for j in vars()["cityComList_" + str(i)]:
                status_value = CBD.nodes[count]['status']
                vars()["city_" + str(i)].nodes[j]['status'] = status_value
                count = count + 1
                
        #Output daily status. Stoped S counting.
        Daily_Result_sum = np.array([0,0,0,0,0,0,0,0]) # S, E, I, R, A, H, Rt, Tt
        
        for i in range(4):
            city_for_func = vars()["city_" + str(i)] 
            Daily_Result_city = Count_status(day,city_for_func)
            Daily_Result_sum = Daily_Result_sum + Daily_Result_city    
            
            print('city_', i, Daily_Result_city)
            
            vars()["cityresult_" + str(i)] = Daily_Result_city
            
        print('day:', day, 'Total', Daily_Result_sum)
        
        Daily_Result_ALL = Daily_Result_sum
        
    elif day > 0 and day < (7*int(day_n/7)):
             
        #time1_str = datetime.datetime.now()
        #print(time1_str)
        
        CBDs=pic.loads(pic.dumps(CBD))
        city_0=pic.loads(pic.dumps(vars()["city_" + str(0)]))
        city_1=pic.loads(pic.dumps(vars()["city_" + str(1)]))
        city_2=pic.loads(pic.dumps(vars()["city_" + str(2)]))
        city_3=pic.loads(pic.dumps(vars()["city_" + str(3)]))
        cityComList_0=pic.loads(pic.dumps(vars()["cityComList_" + str(0)]))
        cityComList_1=pic.loads(pic.dumps(vars()["cityComList_" + str(1)]))
        cityComList_2=pic.loads(pic.dumps(vars()["cityComList_" + str(2)]))
        cityComList_3=pic.loads(pic.dumps(vars()["cityComList_" + str(3)]))
        
        time2_str = datetime.datetime.now()
        print(time2_str)
        
        beta_t=find_beta_t(day,beta_t_S.iloc[day],CBDs,city_0,
                           city_1,
                           city_2,
                           city_3,
                           cityComList_0,
                           cityComList_1,
                           cityComList_2,
                           cityComList_3)
        
        time3_str = datetime.datetime.now()
        print(time3_str)
        
        for t in range(day,day_n):
            beta_t_S.iloc[t]=beta_t

        print('day:',day)
        print('beta_t:',beta_t)    
        #print('beta_t_S:',beta_t_S)
            
        for i in range(4):  #all 4 cities in tau[0] TimeZone
            city_for_func = vars()["city_" + str(i)] # prepare which city to do.
            SEIRAH_SW.SEIRAH_SW(day,city_for_func,beta_t,tau[0])
            
            # Finish tau[0]:life infection. Prepare for tau[1]:life+CBD
            edge_weight_0(city_for_func)
            
        #Interconnect_City2CBD()  
        count = 0
        for i in range(4):
            for j in vars()["cityComList_" + str(i)]:
                status_value = vars()["city_" + str(i)].nodes[j]['status']
                CBD.nodes[count]['status'] = status_value
                count = count + 1 
        
        # SEIRAH process of CBD in a certain day.
        SEIRAH_SW.SEIRAH_SW(day,CBD,beta_t, tau[1])
        
        # SEIRAH precess of cities in working TimeZone of a certain day.    
        for i in range(4):  #all 4 cities in TimeZone tau[1]
            city_for_func = vars()["city_" + str(i)] # prepare which city to do.
            SEIRAH_SW.SEIRAH_SW(day,city_for_func,beta_t,tau[1])
            edge_weight_1(city_for_func)
            
        #Interconnect_CBD2City()      
        count = 0
        for i in range(4):
            for j in vars()["cityComList_" + str(i)]:
                status_value = CBD.nodes[count]['status']
                vars()["city_" + str(i)].nodes[j]['status'] = status_value
                count = count + 1
                
        #Output daily status. Stoped S counting.
        Daily_Result_sum = np.array([0,0,0,0,0,0,0,0]) # S, E, I, R, A, H, Rt, Tt
        
        for i in range(4):
            city_for_func = vars()["city_" + str(i)] 
            Daily_Result_city = Count_status(day,city_for_func)
            Daily_Result_sum = Daily_Result_sum + Daily_Result_city    
            
            print('city_', i, Daily_Result_city)
            
            vars()["cityresult_" + str(i)] = np.vstack((vars()["cityresult_" + str(i)], Daily_Result_city))
            
        print('day:', day, 'Total', Daily_Result_sum)
        
        Daily_Result_ALL = np.vstack((Daily_Result_ALL, Daily_Result_sum))

    else:
        beta_t=beta_t_S.iloc[day]

        for i in range(4):  #all 4 cities in tau[0] TimeZone
            city_for_func = vars()["city_" + str(i)] # prepare which city to do.
            SEIRAH_SW.SEIRAH_SW(day,city_for_func,beta_t,tau[0])
            
            # Finish tau[0]:life infection. Prepare for tau[1]:life+CBD
            edge_weight_0(city_for_func)
            
        #Interconnect_City2CBD()  
        count = 0
        for i in range(4):
            for j in vars()["cityComList_" + str(i)]:
                status_value = vars()["city_" + str(i)].nodes[j]['status']
                CBD.nodes[count]['status'] = status_value
                count = count + 1 
        
        # SEIRAH process of CBD in a certain day.
        SEIRAH_SW.SEIRAH_SW(day,CBD,beta_t, tau[1])
        
        # SEIRAH precess of cities in working TimeZone of a certain day.    
        for i in range(4):  #all 4 cities in TimeZone tau[1]
            city_for_func = vars()["city_" + str(i)] # prepare which city to do.
            SEIRAH_SW.SEIRAH_SW(day,city_for_func,beta_t,tau[1])
            edge_weight_1(city_for_func)
            
        #Interconnect_CBD2City()      
        count = 0
        for i in range(4):
            for j in vars()["cityComList_" + str(i)]:
                status_value = CBD.nodes[count]['status']
                vars()["city_" + str(i)].nodes[j]['status'] = status_value
                count = count + 1
                
        #Output daily status. Stoped S counting.
        Daily_Result_sum = np.array([0,0,0,0,0,0,0,0]) # S, E, I, R, A, H, Rt, Tt
        
        for i in range(4):
            city_for_func = vars()["city_" + str(i)] 
            Daily_Result_city = Count_status(day,city_for_func)
            Daily_Result_sum = Daily_Result_sum + Daily_Result_city    
            
            print('city_', i, Daily_Result_city)
            
            vars()["cityresult_" + str(i)] = np.vstack((vars()["cityresult_" + str(i)], Daily_Result_city))
            
        print('day:', day, 'Total', Daily_Result_sum)
        
        Daily_Result_ALL = np.vstack((Daily_Result_ALL, Daily_Result_sum))

    #print()

Daily_Result_df = pd.DataFrame(Daily_Result_ALL)


#Write GraphML
nx.write_graphml(CBD, "CBD.graphml")
nx.write_graphml(city_0, "City_0.graphml")
nx.write_graphml(city_1, "City_1.graphml")
nx.write_graphml(city_2, "City_2.graphml")
nx.write_graphml(city_3, "City_3.graphml")

#print(cityresult_0)
print(Daily_Result_df)
print('beta_t_S:',beta_t_S)

df = pd.read_csv('new_cases_cr2020.csv', nrows=day_n)
    
real_H = np.zeros((day_n,8))
    
real_H[:,5]=df['H_Shuto']

date_list = df['Date']

day_list = df['Day']

real_H_0 = df['Tokyo']
real_H_1 = df['Kanagawa']
real_H_2 = df['Chiba']
real_H_3 = df['Saitama']
real_H_4 = df['H_Shuto']

#output .csv
#city_0
f_0 = open('output_city_0.csv', "w", newline='')

csv_writer_0 = cs.writer(f_0)

csv_writer_0.writerow(["Date","Day","city_0_beta_t","city_0_S","city_0_E",
                     "city_0_I","city_0_R","city_0_A","city_0_H",
                     "city_0_Rt_NWK","city_0_Tt_NWK","city_0_realH",
                     "city_0_N","city_0_k","city_0_p","city_0_NCom"])

#S = N-E-I-R-A-H
for k in range(day_n):
    array_0 = (date_list[k],day_list[k],beta_t_S[k],
               (para_0['N']-cityresult_0[k,1]-cityresult_0[k,2]-cityresult_0[k,3]-cityresult_0[k,4]-cityresult_0[k,5]),
               cityresult_0[k,1],
             cityresult_0[k,2],cityresult_0[k,3],cityresult_0[k,4],cityresult_0[k,5],
             cityresult_0[k,6],cityresult_0[k,7],real_H_0[k],para_0['N'],
             para_0['k'],para_0['p'],para_0['NCom'])
    csv_writer_0.writerow(array_0)

f_0.close()

#city_1
f_1 = open('output_city_1.csv', "w", newline='')

csv_writer_1 = cs.writer(f_1)

csv_writer_1.writerow(["Date","Day","city_1_beta_t","city_1_S","city_1_E",
                     "city_1_I","city_1_R","city_1_A","city_1_H",
                     "city_1_Rt_NWK","city_1_Tt_NWK","city_1_realH",
                     "city_1_N","city_1_k","city_1_p","city_1_NCom"])

for k in range(day_n):
    array_1 = (date_list[k],day_list[k],beta_t_S[k],
               (para_1['N']-cityresult_1[k,1]-cityresult_1[k,2]-cityresult_1[k,3]-cityresult_1[k,4]-cityresult_1[k,5]),
               cityresult_1[k,1],
             cityresult_1[k,2],cityresult_1[k,3],cityresult_1[k,4],cityresult_1[k,5],
             cityresult_1[k,6],cityresult_1[k,7],real_H_1[k],para_1['N'],
             para_1['k'],para_1['p'],para_1['NCom'])
    csv_writer_1.writerow(array_1)

f_1.close()

#city_2
f_2 = open('output_city_2.csv', "w", newline='')

csv_writer_2 = cs.writer(f_2)

csv_writer_2.writerow(["Date","Day","city_2_beta_t","city_2_S","city_2_E",
                     "city_2_I","city_2_R","city_2_A","city_2_H",
                     "city_2_Rt_NWK","city_2_Tt_NWK","city_2_realH",
                     "city_2_N","city_2_k","city_2_p","city_2_NCom"])

for k in range(day_n):
    array_2 = (date_list[k],day_list[k],beta_t_S[k],
               (para_0['N']-cityresult_2[k,1]-cityresult_2[k,2]-cityresult_2[k,3]-cityresult_2[k,4]-cityresult_2[k,5]),
               cityresult_2[k,1],
             cityresult_2[k,2],cityresult_2[k,3],cityresult_2[k,4],cityresult_2[k,5],
             cityresult_2[k,6],cityresult_2[k,7],real_H_2[k],para_2['N'],
             para_2['k'],para_2['p'],para_2['NCom'])
    csv_writer_2.writerow(array_2)

f_2.close()

#city_3
f_3 = open('output_city_3.csv', "w", newline='')

csv_writer_3 = cs.writer(f_3)

csv_writer_3.writerow(["Date","Day","city_3_beta_t","city_3_S","city_3_E",
                     "city_3_I","city_3_R","city_3_A","city_3_H",
                     "city_3_Rt_NWK","city_3_Tt_NWK","city_3_realH",
                     "city_3_N","city_3_k","city_3_p","city_3_NCom"])

for k in range(day_n):
    array_3 = (date_list[k],day_list[k],beta_t_S[k],
               (para_3['N']-cityresult_3[k,1]-cityresult_3[k,2]-cityresult_3[k,3]-cityresult_3[k,4]-cityresult_3[k,5]),
               cityresult_3[k,1],
             cityresult_3[k,2],cityresult_3[k,3],cityresult_3[k,4],cityresult_3[k,5],
             cityresult_3[k,6],cityresult_3[k,7],real_H_3[k],para_3['N'],
             para_3['k'],para_3['p'],para_3['NCom'])
    csv_writer_3.writerow(array_3)

f_3.close()


#city_CBD
f_4 = open('output_Shuto.csv', "w", newline='')

csv_writer_4 = cs.writer(f_4)

csv_writer_4.writerow(["Date","Day","city_Shuto_beta_t","city_Shuto_S","city_Shuto_E",
                     "city_Shuto_I","city_Shuto_R","city_Shuto_A","city_Shuto_H",
                     "city_Shuto_Rt_NWK","city_Shuto_Tt_NWK","city_Shuto_realH",
                     "city_Shuto_N","city_Shuto_k","city_Shuto_p","city_Shuto_NCom"])

for k in range(day_n):
    array_4 = (date_list[k],day_list[k],beta_t_S[k],
               ((para_0['N']+para_1['N']+para_2['N']+para_3['N'])-Daily_Result_ALL[k,1]-Daily_Result_ALL[k,2]-Daily_Result_ALL[k,3]-Daily_Result_ALL[k,4]-Daily_Result_ALL[k,5]),
               Daily_Result_ALL[k,1],
             Daily_Result_ALL[k,2],Daily_Result_ALL[k,3],Daily_Result_ALL[k,4],
             Daily_Result_ALL[k,5],
             Daily_Result_ALL[k,6],Daily_Result_ALL[k,7],real_H_4[k],(para_0['N']+para_1['N']+para_2['N']+para_3['N']),"NaN","NaN",(para_0['NCom']+para_1['NCom']+para_2['NCom']+para_3['NCom']))
    csv_writer_4.writerow(array_4)

f_4.close()

#S, E, I, R, A, H, Rt, Tt

datem=Daily_Result_df

datem.index=pd.to_datetime(df['Date'].head(day_n))

plt1.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b,%d'))
#plt1.gca().xaxis.set_major_locator(mdates.WeekdayLocator()) #
plt1.gca().xaxis.set_major_locator(mdates.MonthLocator()) #

#plt.plot(Daily_Result_ALL[:,0],color ='#83C69F',label = 'Susceptible')
plt1.plot(datem[1],color ='#F8D200',label = 'Exposed')
plt1.plot(datem[2],color ='#ED6E68',label = 'Infected')
#plt.plot(Daily_Result_ALL[:,3],color ='#009DAA',label = 'Recovered')
plt1.plot(datem[4],color ='#6A9FD3',label = 'Asymptomatic')
plt1.plot(datem[5],color ='#7A1E71',label = 'Hospitalized')

plt1.title('') #Interconnected SEIRAH_SW Model
plt1.legend()
#plt1.xlabel('Day')

plt1.gcf().autofmt_xdate()

plt1.ylabel('Number')
#plt1.figure()
plt1.savefig('Model')
plt1.figure()

plt2.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b,%d'))
#plt2.gca().xaxis.set_major_locator(mdates.WeekdayLocator()) #
plt2.gca().xaxis.set_major_locator(mdates.MonthLocator()) #

dateb=beta_t_S

dateb.index=pd.to_datetime(df['Date'].head(day_n))

plt2.plot(dateb,color ='#7A2871',label = '') #beta_t

plt2.gcf().autofmt_xdate()

plt2.title('') #Beta_t
plt2.legend()
#plt2.xlabel('Day')
plt2.ylabel('Number')
plt2.savefig('Beta_t')
plt2.figure()

plt3.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b,%d'))
#plt3.gca().xaxis.set_major_locator(mdates.WeekdayLocator()) #
plt3.gca().xaxis.set_major_locator(mdates.MonthLocator()) #

dateH=pd.DataFrame(real_H[:,5])

dateH.index=pd.to_datetime(df['Date'].head(day_n))

plt3.plot(datem[5],color ='#7A1E71',label = 'Predicted Hospitalized')
plt3.plot(dateH,color ='#6A9FD3',label = 'Real Hospitalized')

plt3.gcf().autofmt_xdate()

plt3.title('') #Hospitalized
plt3.legend()
#plt3.xlabel('Day')
plt3.ylabel('Number')
plt3.savefig('Hospitalized')
plt3.figure()
