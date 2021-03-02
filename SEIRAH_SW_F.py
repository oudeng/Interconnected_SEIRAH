#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 01:29:34 2021
@author: Deng, Sandy

Count status of SEIRAH in target network. 

Note:
    Daily data. NOT for DayHours.
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

def SEIRAH_SW_F(_nwks,beta_t,tau):
    
    para_Epidemic = {'sigma':0.2,'p1':0.18,'p2':0.3,'l_AH':0.05,'l_IH':0.3,
                 'g_AR':0.07,'g_HR':0.1}
    
    _nwks_F=pic.loads(pic.dumps(_nwks))

    for i in range(_nwks_F.number_of_nodes()):
        
        # Pick node status. Infect neighbours by probability.        
        if _nwks_F.nodes[i]['status'] == 'infe':
            # Precisely, infection first or change first probabilistically 0.5:0.5         
            if random.random() > 0.5:  #Infecte firstly pattern.            
                # I infect S to E probabilistically
                threshold = 1 - beta_t*tau
                for nbr in _nwks_F[i]:
                    probability = random.random()*_nwks_F.edges[i,nbr]['weight']
                    if _nwks_F.nodes[nbr]['status'] == 'susc' and probability > threshold:
                        _nwks_F.nodes[nbr]['status'] = 'expo'
                        #_nwks.nodes[nbr]['E_1stday'] = time_stamp
                        _nwks_F.nodes[i]['Infe_other'] = _nwks_F.nodes[i]['Infe_other'] +1
                    #print(_nwks.nodes[nbr],_nwks.nodes[i])
                    
                # I changes to H
                threshold = 1 - para_Epidemic['l_IH']*tau
                if random.random() > threshold:
                    _nwks_F.nodes[i]['status'] = 'hosp'
                    #_nwks.nodes[i]['H_1stday'] = time_stamp
            
            else:   #while random.random()<0.5. Change first pattern.             
                # I changes to H
                threshold = 1 - para_Epidemic['l_IH']*tau
                if random.random() > threshold:
                    _nwks_F.nodes[i]['status'] = 'hosp'
                    #_nwks.nodes[i]['H_1stday'] = time_stamp
                    
                else:
                    # if no change, I infect S to E probabilistically
                    threshold = 1 - beta_t*tau
                    for nbr in _nwks_F[i]:
                        probability = random.random()*_nwks_F.edges[i,nbr]['weight']
                        if _nwks_F.nodes[nbr]['status'] == 'susc' and probability > threshold:
                            _nwks_F.nodes[nbr]['status'] = 'expo'
                            #_nwks.nodes[nbr]['E_1stday'] = time_stamp
                            _nwks_F.nodes[i]['Infe_other'] = _nwks_F.nodes[i]['Infe_other'] +1
            #continue  #Continue to next node status process.
            
        if _nwks_F.nodes[i]['status'] == 'asym':          
            # Precisely, infection first or change first probabilistically 0.5:0.5         
            if 0 < random.random() <= 0.33: # Infection first pattern.
                # A infect S to E probabilistically
                threshold = 1- beta_t*tau
                for nbr in _nwks_F[i]:
                    probability = random.random()*_nwks_F.edges[i,nbr]['weight']
                    if _nwks_F.nodes[nbr]['status'] == 'susc' and probability > threshold:
                        _nwks_F.nodes[nbr]['status'] = 'expo'
                        #_nwks.nodes[nbr]['E_1stday'] = time_stamp
                        _nwks_F.nodes[i]['Infe_other'] = _nwks_F.nodes[i]['Infe_other'] +1
        
                # A changes to R
                threshold = 1 - (1-para_Epidemic['p2'])*para_Epidemic['g_AR']*tau
                if random.random() > threshold:
                    _nwks_F.nodes[i]['status'] = 'reco'
                    #_nwks.nodes[i]['R_1stday'] = time_stamp
                    
                else:
                    # A changes to H
                    threshold = 1 - para_Epidemic['p2']*para_Epidemic['l_AH']*tau
                    if random.random() > threshold:
                        _nwks_F.nodes[i]['status'] = 'hosp'
                        #_nwks.nodes[i]['H_1stday'] = time_stamp

            elif 0.33 < random.random() <= 0.67: # Infection first pattern.
                # A infect S to E probabilistically
                threshold = 1- beta_t*tau
                for nbr in _nwks_F[i]:
                    probability = random.random()*_nwks_F.edges[i,nbr]['weight']
                    if _nwks_F.nodes[nbr]['status'] == 'susc' and probability > threshold:
                        _nwks_F.nodes[nbr]['status'] = 'expo'
                        #_nwks.nodes[nbr]['E_1stday'] = time_stamp
                        _nwks_F.nodes[i]['Infe_other'] = _nwks_F.nodes[i]['Infe_other'] +1     
        
                # A changes to H
                threshold = 1 - para_Epidemic['p2']*para_Epidemic['l_AH']*tau
                if random.random() > threshold:
                    _nwks_F.nodes[i]['status'] = 'hosp'
                    #_nwks.nodes[i]['H_1stday'] = time_stamp
                    
                else:
                    # A changes to R
                    threshold = 1 - (1-para_Epidemic['p2'])*para_Epidemic['g_AR']*tau
                    if random.random() > threshold:
                        _nwks_F.nodes[i]['status'] = 'reco'
                        #_nwks.nodes[i]['R_1stday'] = time_stamp
                        
            elif 0.67 < random.random() <= 1.0:      
                # A changes to R
                threshold = 1 - (1-para_Epidemic['p2'])*para_Epidemic['g_AR']*tau
                if random.random() > threshold:
                    _nwks_F.nodes[i]['status'] = 'reco'
                    #_nwks.nodes[i]['R_1stday'] = time_stamp
                    continue
                    
                else:
                    # A changes to H
                    threshold = 1 - para_Epidemic['p2']*para_Epidemic['l_AH']*tau
                    if random.random() > threshold:
                        _nwks_F.nodes[i]['status'] = 'hosp'
                        #_nwks.nodes[i]['H_1stday'] = time_stamp
             
                # A infect S to E probabilistically
                if _nwks_F.nodes[i]['status'] == 'asym':
                    threshold = 1- beta_t*tau
                    for nbr in _nwks_F[i]:
                        probability = random.random()*_nwks_F.edges[i,nbr]['weight']
                        if _nwks_F.nodes[nbr]['status'] == 'susc' and probability > threshold:
                            _nwks_F.nodes[nbr]['status'] = 'expo'
                           # _nwks.nodes[nbr]['E_1stday'] = time_stamp
                            _nwks_F.nodes[i]['Infe_other'] = _nwks_F.nodes[i]['Infe_other'] +1

            #continue
                
        if _nwks_F.nodes[i]['status'] == 'expo':
            # Precisely, infection first or change first probabilistically 0.5:0.5         
            if random.random() > 0.5:                

                # E infect S to E probabilistically
                threshold = 1 - beta_t*tau
                for nbr in _nwks_F[i]:
                    probability = random.random()*_nwks_F.edges[i,nbr]['weight']
                    if _nwks_F.nodes[nbr]['status'] == 'susc' and probability > threshold:
                        _nwks_F.nodes[nbr]['status'] = 'expo'
                        #_nwks.nodes[nbr]['E_1stday'] = time_stamp
                        _nwks_F.nodes[i]['Infe_other'] = _nwks_F.nodes[i]['Infe_other'] +1
                               
                # E changes to I
                threshold = 1 - (1-para_Epidemic['p1'])*para_Epidemic['sigma']*tau
                if random.random() > threshold:
                    _nwks_F.nodes[i]['status'] = 'infe'
                    #_nwks.nodes[i]['I_1stday'] = time_stamp
                    
                else:
                    # E changes to A
                    threshold = 1 - para_Epidemic['p1']*para_Epidemic['sigma']*tau
                    if random.random() > threshold:
                        _nwks_F.nodes[i]['status'] = 'asym'
                        #_nwks.nodes[i]['A_1stday'] = time_stamp               
            else:
                # E changes to I
                threshold = 1 - (1-para_Epidemic['p1'])*para_Epidemic['sigma']*tau
                if random.random() > threshold:
                    _nwks_F.nodes[i]['status'] = 'infe'
                    #_nwks.nodes[i]['I_1stday'] = time_stamp
                    
                else:
                    # E changes to A
                    threshold = 1 - para_Epidemic['p1']*para_Epidemic['sigma']*tau
                    if random.random() > threshold:
                        _nwks_F.nodes[i]['status'] = 'asym'
                        #_nwks.nodes[i]['A_1stday'] = time_stamp
             
                # E infect S to E probabilistically
                threshold = 1 - beta_t*tau
                for nbr in _nwks_F[i]:
                    probability = random.random()*_nwks_F.edges[i,nbr]['weight']
                    if _nwks_F.nodes[nbr]['status'] == 'susc' and probability > threshold:
                        _nwks_F.nodes[nbr]['status'] = 'expo'
                        #_nwks.nodes[nbr]['E_1stday'] = time_stamp
                        _nwks_F.nodes[i]['Infe_other'] = _nwks_F.nodes[i]['Infe_other'] +1
            #continue  
                                   
        if _nwks_F.nodes[i]['status'] == 'hosp':
            # H changes to R
            threshold = 1 - para_Epidemic['g_HR']*tau
            if random.random() > threshold:
                _nwks_F.nodes[i]['status'] = 'reco'
                #_nwks.nodes[i]['R_1stday'] = time_stamp   
            #continue

    return _nwks_F