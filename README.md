# Intconnected_SEIRAH
An Extended Epidemic Model on Interconnected Networks for COVID-19 to Explore the Epidemic Dynamics

Purpose:
    Experiment attached to our research paper titled “An Extended Epidemic Model
    on Interconnected Networks for COVID-19 to Explore the Epidemic Dynamics”. 
    This paper was coauthored by Ou Deng, Kiichi Tago, and Qun Jin. 
    Networked Information Systems Laboratory (NISL)
    Faculty of Human Sciences, Waseda University

Environment:

    1. Software: Python 3.7 with NetworkX library.
    2. Hardware: Alicloud (high performance computing is necessary)
    
Experiment package includes:

    1. SEIRAH_main.py: main program.
    2. SEIRAH_PredictEnding.py: predict COVID-19 ending date.
    3. SEIRAH_SW_F.py: function imported.
    4. SEIRAH_SW.py: function imported.
    5. new_cases_cr2020.csv: Dataset.

Dataset:

    1. new_cases_cr2020.csv. More details mentioned in the above paper.
        - Daily new cases time series of Tokyo, Kanagawa, Saitama, and Chiba,
          since COVID-19 outbreak on Feb., 2020 in Japan.
        - Daily commuting ratio time series, since Mar.,2020 Japan officialy
          started to announce it in public.
    2. Fixed epidemiological parameter setting according to COVID-19 situation in 2020.
       It is written in python for this version.

What 'Interconnected SEIRAH' can do:

    1. Realize Algorithm 1 and 2 of the above mentioned paper.
    2. Simulate SEIRAH transition in time series (csv output).
    3. Verify netwrok topology by adjusting parameter setting, especially p.
    4. Predict propagation ending date.

Limitations of 'Interconnected SEIRAH':

    1. Target social networks. It is suitable for metropolitan scale simulation. If much more extensive, such as national scale, the intervention estimator results become meaningless. If much smaller, the stochastic processing results will not reliable enough.
    2. Computating time. Inside causal inference processes cost much time.
    
Its Milestones:

    Oct, 29, 2020  Full functions combined with cloud computing outputs.
    Nov, 16, 2020  360k-node network varietions verification.
    Feb, 2,  2021  Final combination of 'main-module' Python pacakage.
    
Contact person:

    Ou DENG
    Email: dengou@toki.waseda.jp
    Networked Information Systems Laboratory
    (NISLAB: https://nislab.human.waseda.ac.jp/)
    Graduate School of Human Sciences, Waseda University
    早稲田大学院　人間科学研究科
