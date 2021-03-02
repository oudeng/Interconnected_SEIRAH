#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 01:29:34 2021
@author: Deng, Sandy

Predict the infection ending:
    SEIRAH time series proceeding from the last situation of SEIRAH_main.py,
    until the condition of E+A+I=0.
"""

import networkx as nx
import random
from random import sample
import matplotlib.pyplot as plt1
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import copy as co
import csv as cs
import datetime
import SEIRAH_SW as SEIRAH_SW

random.seed(2020)

"""
Network Generation and Initialization.
NetworkX official documents:
https://networkx.github.io/documentation/stable/_modules/networkx/generators/random_graphs.html#newman_watts_strogatz_graph
City networks generated. H be isolated(node weight=0). 1stday of each status.
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
    # Record *_1stday. S no need to record.
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
Count statuses of SEIRAH in target network. Daily data. NOT for TimeZone.
"""

def Count_status(__time_stamp, __nwk):
    
    num_S = num_E = num_I = num_R = num_A = num_H = 0
    count_Infe_other = count_T = count_H = 0
    Rt = Tt = 0

    for i in range(__nwk.number_of_nodes()):
#        S is big number. No count to save computing time.
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
            count_H = count_H +1
            count_Infe_other = count_Infe_other + __nwk.nodes[i]['Infe_other']
            #Count how many days from E to H
            count_T = count_T + __nwk.nodes[i]['H_1stday'] - __nwk.nodes[i]['E_1stday']

    if count_H ==0 or num_H == 0:  #Avoid error while H=0
        Rt = Tt = 0
    else:
        Rt = count_Infe_other / count_H
        Tt = count_T / count_H
        #Tt = format(Rt / beta_t,'.2f')
    
    __Count_status = [num_S, num_E, num_I, num_R, num_A, num_H, Rt, Tt]
    
    return __Count_status

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #       Above sections are FUNCTIONS.             # #
# #       MAIN() is as follows.                     # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

para_Epidemic = {'sigma':0.2,'p1':0.18,'p2':0.3,'l_AH':0.05,'l_IH':0.3,
                 'g_AR':0.07,'g_HR':0.1}

#beta = [0.1,0.1,0.1] # B_IS＝0.3, B_AS＝0.2, B_ES＝0.2 as default
#beta_t = 0.1 # by Sandy:

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

#Import parameters from .csv, by Sandy

df0 = pd.read_csv('output_city_0.csv')
df1 = pd.read_csv('output_city_1.csv')
df2 = pd.read_csv('output_city_2.csv')
df3 = pd.read_csv('output_city_3.csv')
df_Shuto = pd.read_csv('output_Shuto.csv')

last_date = df_Shuto['Date'].values[-1]
predict_start_date=pd.Series([pd.to_datetime(last_date)])

df_RealH = pd.read_csv('new_cases_cr2020.csv').set_index('Date')
Real_H_Shuto = df_RealH.at[last_date,'H_Shuto']

Real_H_0 = df_RealH.at[last_date,'Tokyo']
Real_H_1 = df_RealH.at[last_date,'Kanagawa']
Real_H_2 = df_RealH.at[last_date,'Chiba']
Real_H_3 = df_RealH.at[last_date,'Saitama']

H_Ratio_0 = Real_H_0 / Real_H_Shuto
H_Ratio_1 = Real_H_1 / Real_H_Shuto
H_Ratio_2 = Real_H_2 / Real_H_Shuto
H_Ratio_3 = Real_H_3 / Real_H_Shuto

#City_0
#average beta
ave=6
beta_sum = 0
while ave >= 0:
    beta_sum = beta_sum + df_Shuto['city_Shuto_beta_t'].values[-1-ave]
    print(beta_sum)
    ave -=1
    
    
beta = beta_sum/7 # Prediction period. 7(days) is as parameter.

print("Start_Beta:",beta)  #Original avg+1

CR = 0.75 # Predict

print(predict_start_date)

#expo+infe+asym+hosp+reco

para_0 = {'N':int(df0['city_0_N'].values[-1]),
          'k':int(df0['city_0_k'].values[-1]),
          'p':df0['city_0_p'].values[-1],
          'NCom':df0['city_0_NCom'].values[-1],
          'expo':int(df_Shuto['city_Shuto_E'].values[-1]*H_Ratio_0),
          'infe':int(df_Shuto['city_Shuto_I'].values[-1]*H_Ratio_0),
          'asym':int(df_Shuto['city_Shuto_A'].values[-1]*H_Ratio_0),
          'hosp':int(df_Shuto['city_Shuto_H'].values[-1]*H_Ratio_0),
          'reco':int(df_Shuto['city_Shuto_R'].values[-1]*H_Ratio_0)}

para_1 = {'N':int(df1['city_1_N'].values[-1]),
          'k':int(df1['city_1_k'].values[-1]),
          'p':df1['city_1_p'].values[-1],
          'NCom':df1['city_1_NCom'].values[-1],
          'expo':int(df_Shuto['city_Shuto_E'].values[-1]*H_Ratio_1),
          'infe':int(df_Shuto['city_Shuto_I'].values[-1]*H_Ratio_1),
          'asym':int(df_Shuto['city_Shuto_A'].values[-1]*H_Ratio_1),
          'hosp':int(df_Shuto['city_Shuto_H'].values[-1]*H_Ratio_1),
          'reco':int(df_Shuto['city_Shuto_R'].values[-1]*H_Ratio_1)}

para_2 = {'N':int(df2['city_2_N'].values[-1]),
          'k':int(df2['city_2_k'].values[-1]),
          'p':df2['city_2_p'].values[-1],
          'NCom':df2['city_2_NCom'].values[-1],
          'expo':int(df_Shuto['city_Shuto_E'].values[-1]*H_Ratio_2),
          'infe':int(df_Shuto['city_Shuto_I'].values[-1]*H_Ratio_2),
          'asym':int(df_Shuto['city_Shuto_A'].values[-1]*H_Ratio_2),
          'hosp':int(df_Shuto['city_Shuto_H'].values[-1]*H_Ratio_2),
          'reco':int(df_Shuto['city_Shuto_R'].values[-1]*H_Ratio_2)}

para_3 = {'N':int(df3['city_3_N'].values[-1]),
          'k':int(df3['city_3_k'].values[-1]),
          'p':df3['city_3_p'].values[-1],
          'NCom':df3['city_3_NCom'].values[-1],
          'expo':int(df_Shuto['city_Shuto_E'].values[-1]*H_Ratio_3),
          'infe':int(df_Shuto['city_Shuto_I'].values[-1]*H_Ratio_3),
          'asym':int(df_Shuto['city_Shuto_A'].values[-1]*H_Ratio_3),
          'hosp':int(df_Shuto['city_Shuto_H'].values[-1]*H_Ratio_3),
          'reco':int(df_Shuto['city_Shuto_R'].values[-1]*H_Ratio_3)}

print(para_0,para_1,para_2,para_3)

para_N_0=co.copy(para_0)
para_N_1=co.copy(para_1)
para_N_2=co.copy(para_2)
para_N_3=co.copy(para_3)

CBD_N = para_0['NCom']+para_1['NCom']+para_2['NCom']+para_3['NCom']
CBD_k = 8   
CBD_p = 0.05

cityComList =[]
l=0

Daily_Result_sum = np.array([df_Shuto['city_Shuto_S'].values[-1],df_Shuto['city_Shuto_E'].values[-1],
                             df_Shuto['city_Shuto_I'].values[-1],df_Shuto['city_Shuto_R'].values[-1],
                             df_Shuto['city_Shuto_A'].values[-1],df_Shuto['city_Shuto_H'].values[-1],
                             df_Shuto['city_Shuto_Rt_NWK'].values[-1],df_Shuto['city_Shuto_Tt_NWK'].values[-1]]) 
# Initial S,E,I,R,A,H,Rt,Tt

## parameters input END ######

'''
Variable notes:
    nwk: network. 函数形参（parameter），（实参=argument）。
    city_*: citys in simulation. city_0 is CENTER, other are SKIRT(s).
    cityComList_*: commuting nodes index list of city_*.
    cityComList: sequenced combination of cityComList_0,1,2,3.
                for interconnection function.
'''

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


# Especially while periodcally oscilation ocours, change to fixed period.

#Initial Status
day = 0

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

# Initial day, β, by Sandy

print('Start Ending Predict:')

while ((Daily_Result_sum[1]+Daily_Result_sum[2]+Daily_Result_sum[4]) > 0 and day < 1000):
    
    predict_start_date.loc[day + 1]=predict_start_date.iloc[day] + datetime.timedelta(days=1)
    
    day = day + 1
    
    count = 0  # CBD nodes index.
    
    for cr in range(4):  # 4 cities.
    
        para = vars()["para_" + str(cr)]
        para_N = vars()["para_N_" + str(cr)]
        
        #print(para,para_N)

        para['NCom'] = int(para_N['NCom']*CR)
        
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
    
    for i in range(4):  #all 4 cities in tau[0] TimeZone
        city_for_func = vars()["city_" + str(i)] # prepare which city to do.
        SEIRAH_SW.SEIRAH_SW(day,city_for_func,beta,tau[0])
        
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
    SEIRAH_SW.SEIRAH_SW(day,CBD,beta, tau[1])
    
    # SEIRAH precess of cities in working TimeZone of a certain day.    
    for i in range(4):  #all 4 cities in TimeZone tau[1]
        city_for_func = vars()["city_" + str(i)] # prepare which city to do.
        SEIRAH_SW.SEIRAH_SW(day,city_for_func,beta,tau[1])
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
#End SEIRAH, By Sandy

#Write GraphML
nx.write_graphml(CBD, "CBD.graphml")
nx.write_graphml(vars()["city_" + str(0)], "City_0_Predict.graphml")
nx.write_graphml(vars()["city_" + str(1)] , "City_1_Predict.graphml")
nx.write_graphml(vars()["city_" + str(2)] , "City_2_Predict.graphml")
nx.write_graphml(vars()["city_" + str(3)] , "City_3_Predict.graphml")

#print(cityresult_0)
print(Daily_Result_df)
print('beta:',beta)
print('day:',day)

date_list = predict_start_date

day_list = '0'

real_H_list = '0'

#output .csv
#city_0
f_0 = open('predict_output_city_0.csv', "w", newline='')

csv_writer_0 = cs.writer(f_0)

csv_writer_0.writerow(["Date","Day","city_0_beta_t","city_0_S","city_0_E",
                     "city_0_I","city_0_R","city_0_A","city_0_H",
                     "city_0_Rt_NWK","city_0_Tt_NWK","city_0_realH",
                     "city_0_N","city_0_k","city_0_p","city_0_Ncom"])

#S = N-E-I-R-A-H
for k in range(day):
    array_0 = (date_list[k],day,beta,
               (para_0['N']-cityresult_0[k,1]-cityresult_0[k,2]-cityresult_0[k,3]-cityresult_0[k,4]-cityresult_0[k,5]),
               cityresult_0[k,1],
             cityresult_0[k,2],cityresult_0[k,3],cityresult_0[k,4],cityresult_0[k,5],
             cityresult_0[k,6],cityresult_0[k,7],"NaN",para_0['N'],
             para_0['k'],para_0['p'],para_0['NCom'])
    csv_writer_0.writerow(array_0)

f_0.close()

#city_1
f_1 = open('predict_output_city_1.csv', "w", newline='')

csv_writer_1 = cs.writer(f_1)

csv_writer_1.writerow(["Date","Day","city_1_beta_t","city_1_S","city_1_E",
                     "city_1_I","city_1_R","city_1_A","city_1_H",
                     "city_1_Rt_NWK","city_1_Tt_NWK","city_1_realH",
                     "city_1_N","city_1_k","city_1_p","city_1_Ncom"])

for k in range(day):
    array_1 = (date_list[k],day,beta,
               (para_1['N']-cityresult_1[k,1]-cityresult_1[k,2]-cityresult_1[k,3]-cityresult_1[k,4]-cityresult_1[k,5]),
               cityresult_1[k,1],
             cityresult_1[k,2],cityresult_1[k,3],cityresult_1[k,4],cityresult_1[k,5],
             cityresult_1[k,6],cityresult_1[k,7],"NaN",para_1['N'],
             para_1['k'],para_1['p'],para_1['NCom'])
    csv_writer_1.writerow(array_1)

f_1.close()

#city_2
f_2 = open('predict_output_city_2.csv', "w", newline='')

csv_writer_2 = cs.writer(f_2)

csv_writer_2.writerow(["Date","Day","city_2_beta_t","city_2_S","city_2_E",
                     "city_2_I","city_2_R","city_2_A","city_2_H",
                     "city_2_Rt_NWK","city_2_Tt_NWK","city_2_realH",
                     "city_2_N","city_2_k","city_2_p","city_2_Ncom"])

for k in range(day):
    array_2 = (date_list[k],day,beta,
               (para_0['N']-cityresult_2[k,1]-cityresult_2[k,2]-cityresult_2[k,3]-cityresult_2[k,4]-cityresult_2[k,5]),
               cityresult_2[k,1],
             cityresult_2[k,2],cityresult_2[k,3],cityresult_2[k,4],cityresult_2[k,5],
             cityresult_2[k,6],cityresult_2[k,7],"NaN",para_2['N'],
             para_2['k'],para_2['p'],para_2['NCom'])
    csv_writer_2.writerow(array_2)

f_2.close()

#city_3
f_3 = open('predict_output_city_3.csv', "w", newline='')

csv_writer_3 = cs.writer(f_3)

csv_writer_3.writerow(["Date","Day","city_3_beta_t","city_3_S","city_3_E",
                     "city_3_I","city_3_R","city_3_A","city_3_H",
                     "city_3_Rt_NWK","city_3_Tt_NWK","city_3_realH",
                     "city_3_N","city_3_k","city_3_p","city_3_Ncom"])

for k in range(day):
    array_3 = (date_list[k],day,beta,
               (para_3['N']-cityresult_3[k,1]-cityresult_3[k,2]-cityresult_3[k,3]-cityresult_3[k,4]-cityresult_3[k,5]),
               cityresult_3[k,1],
             cityresult_3[k,2],cityresult_3[k,3],cityresult_3[k,4],cityresult_3[k,5],
             cityresult_3[k,6],cityresult_3[k,7],"NaN",para_3['N'],
             para_3['k'],para_3['p'],para_3['NCom'])
    csv_writer_3.writerow(array_3)

f_3.close()


#city_CBD
f_4 = open('predict_output_Shuto.csv', "w", newline='')

csv_writer_4 = cs.writer(f_4)

csv_writer_4.writerow(["Date","Day","city_Shuto_beta_t","city_Shuto_S","city_Shuto_E",
                     "city_Shuto_I","city_Shuto_R","city_Shuto_A","city_Shuto_H",
                     "city_Shuto_Rt_NWK","city_Shuto_Tt_NWK","city_Shuto_realH",
                     "city_Shuto_N","city_Shuto_k","city_Shuto_p","city_Shuto_Ncom"])

for k in range(day):
    array_4 = (date_list[k],day,beta,
               ((para_0['N']+para_1['N']+para_2['N']+para_3['N'])-Daily_Result_ALL[k,1]-Daily_Result_ALL[k,2]-Daily_Result_ALL[k,3]-Daily_Result_ALL[k,4]-Daily_Result_ALL[k,5]),
               Daily_Result_ALL[k,1],
             Daily_Result_ALL[k,2],Daily_Result_ALL[k,3],Daily_Result_ALL[k,4],
             Daily_Result_ALL[k,5],
             Daily_Result_ALL[k,6],Daily_Result_ALL[k,7],"NaN",(para_0['N']+para_1['N']+para_2['N']+para_3['N']),"NaN","NaN",(para_0['NCom']+para_1['NCom']+para_2['NCom']+para_3['NCom']))
    csv_writer_4.writerow(array_4)

f_4.close()

# Interconnected SEIRAH_SW. Considering S and R are bigger than others, 
# Comment-out S and R to show E-I-A-H more clearly.

#S, E, I, R, A, H, Rt, Tt

#print('Test:',df['Date'].head(day_n))

#print(datem[1])

datem=Daily_Result_df

datem.index=predict_start_date

plt1.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b,%d'))
#plt1.gca().xaxis.set_major_locator(mdates.WeekdayLocator()) #
plt1.gca().xaxis.set_major_locator(mdates.MonthLocator()) #

#plt.plot(Daily_Result_ALL[:,0],color ='#83C69F',label = 'Susceptible')
plt1.plot(datem[1],color ='#F8D200',label = 'Exposed')
plt1.plot(datem[2],color ='#ED6E68',label = 'Infected')
#plt.plot(Daily_Result_ALL[:,3],color ='#009DAA',label = 'Recovered')
plt1.plot(datem[4],color ='#6A9FD3',label = 'Asymptomatic')
plt1.plot(datem[5],color ='#7A1E71',label = 'Hospitalized')

plt1.title('Interconnected SEIRAH_SW Model_Predict') #Interconnected SEIRAH_SW Model
plt1.legend()
#plt1.xlabel('Day')

plt1.gcf().autofmt_xdate()

plt1.ylabel('Number')
#plt1.figure()
plt1.savefig('Model_Predict')
plt1.figure()
