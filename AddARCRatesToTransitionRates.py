"""
a script to determine ratios of fluorescence at different wavelengths as a function of temperature
"
@author: ebn
"""

import numpy as np
import pandas as pd
import os
import math
from scipy.special import wofz
from scipy.constants import k as C_k
from scipy.constants import m_e as m_e
from scipy.constants import m_p as m_p
from scipy.constants import hbar as hbar
from scipy import constants as cts
from fractions import Fraction
import arc

species = 'Rb'  #Currently you may choose 'Rb' or 'Sr'
atom = arc.Rubidium()
new_file='Rb1-transition-rates_ARC_append.xlsx'



os.chdir(r'C:\Users\ebn1\OneDrive - NIST\CoBRAS\Fluorescence-Ratio-Calculator')  #Change this to the directory of the .py file

#edit the location of the xlsx file containing the transition rates for the desired species, assumes these files are also in the current working directory, else do appropriate os operations
transition_file ={
    'Sr':  'Sr1-transition-rates.xlsx',
    'Rb':  'Rb1-transition-rates.xlsx',
    }

state_file ={
    'Sr': 'Sr-states.xlsx',
    'Rb': 'Rb-states.xlsx'}



headers = ['Initial Configuration',
           'Initial term',
           'Initial J',
           'Final configuration',
           'Final term',
           'Final J',
           'Wavelength (nm)',
           'Wavelength uncertainty (nm)',
           'Wavelength w/uncertainty (nm)',
           'Matrix element (a.u)',
           'Matrix el. uncertainty (a.u)',
           'Matrix el. w/uncertainty (a.u)',
           'Transition rate (s-1)',
           'Transition rate uncertainty',
           'Transition rate w/uncertainty (s-1)',
           'Branching ratio',
           'Branching ratio uncertainty',
           'Branching ratio w/uncertainty',
           'Lifetime reference',
           'Matrix Element Reference'
]


dtypes ={'Initial Configuration': 'str',
           'Initial term': 'str',
           'Initial J': 'int',
           'Final configuration': 'str',
           'Final term': 'str',
           'Final J': 'int',
           'Wavelength (nm)': 'float',
           'Wavelength uncertainty (nm)': 'float',
           'Wavelength w/uncertainty (nm)': 'str',
           'Matrix element (a.u)' :'float',
           'Matrix el. uncertainty (a.u)':'float',
           'Matrix el. w/uncertainty (a.u)': 'str',
           'Transition rate (s-1)': 'float',
           'Transition rate uncertainty': 'float',
           'Transition rate w/uncertainty (s-1)' :'str',
           'Branching ratio': 'float',
           'Branching ratio uncertainty': 'float',
           'Branching ratio w/uncertainty' :'str',
           'Lifetime reference':'str',
           'Matrix Element Reference': 'str'
}




transition_data = pd.read_excel(transition_file[species])
state_data = pd.read_excel(state_file[species])

#remove any unnecessary spaces from the 'J' string
if isinstance(state_data['J'].iloc[0], str):
    state_data['J']  = state_data['J'].str.strip()

    transition_data['Initial J']  = transition_data['Initial J'].str.strip()

    transition_data['Final J']  = transition_data['Final J'].str.strip()
#%%
def get_index(config, term, j):
    #initial_row = pd.Series({'Configuration': config,
    #               'Term': term,
    #               'J': j })
    #return states_data.apply(lambda row: row.equals(initial_row), axis=1)
    return state_data[( state_data['Configuration']== config) &\
                (state_data['Term']== term )&\
                (state_data['J']== j) ].index[0]
    
#print(get_index('7p', '2P', '3/2'))

######
# Now we have state_data where every state has a unique index.  And whe have transition data, where info about each transition is organized.
# To make it easier to create population transfer matrixes, we want to assign the state index to the initial and final state in transition_data.
# Here, we find the state index for each initial and final state in the transition_data.  Then we append it to transition_data
# This will let us look up transitions more efficiently
######
td_i_list = [] 
td_f_list =[]
for i in range(len(transition_data)):
    td_i_list.append(get_index( 
        transition_data['Initial Configuration'].iloc[i], 
        transition_data['Initial term'].iloc[i], 
        transition_data['Initial J'].iloc[i]
        ))
    
    td_f_list.append(get_index( 
        transition_data['Final configuration'].iloc[i], 
        transition_data['Final term'].iloc[i], 
        transition_data['Final J'].iloc[i]
        ))
#print(td_i_list)              
#print(td_f_list) 
td_i_df = pd.DataFrame({'Initial Index' : td_i_list})
td_f_df = pd.DataFrame({'Final Index' : td_f_list})

if 'Initial Index' not in transition_data.columns:
    transition_data = pd.concat([transition_data, td_i_df],axis = 1)
if 'Final Index' not in transition_data.columns:
    transition_data = pd.concat([transition_data, td_f_df],axis = 1)







def get_Transition_Rate(i):  #a function to check consistency.  Takes in transition_data index idx, calculated the transition rate from the matrix element and wavelength
    wav = transition_data['Wavelength (nm)'].iloc[i] *1e-9
    mu  = transition_data['Matrix element (a.u)'].iloc[i]*(cts.e*cts.value('Bohr radius'))
    Ji = transition_data['Initial J'].iloc[i]
    if isinstance(Ji, str):
        Ji = float(Fraction(Ji))
    Jf = transition_data['Final J'].iloc[i]
    if isinstance(Jf, str):
        Jf = float(Fraction(Jf))
    Rate_portal = transition_data['Transition rate (s-1)'].iloc[i]
    
    rate = (16*np.pi**3 * mu**2)/(3 * cts.epsilon_0*cts.h * wav**3)/(2*Ji+1)
    return rate, Rate_portal, rate/Rate_portal, Ji, Jf

#print(get_Transition_Rate(1))
"""
def calculate_Transition_Rate(i, j, temperature = 0): #determine the transition rate between states indexed i and j.  Account for BBR
    
    #is i -> j a decay?
    if len(transition_data[( transition_data['Initial Index'] == i) &  (transition_data['Final Index'] == j)]) == 1:
        df = transition_data[( transition_data['Initial Index'] == i) &  (transition_data['Final Index'] == j)]
        omega = 2.0 * np.pi * cts.c/  (df['Wavelength (nm)'].iloc[0]*1e-9)
        gamma =df['Transition rate (s-1)'].iloc[0]
        
        modeOccupationTerm = 1
        degeneracyTerm = 1

        # only possible by absorbing thermal photons ?
        if (hbar * omega < 100 * C_k * temperature) and (omega > 1e2):
            modeOccupationTerm += 1. / \
                (np.exp(hbar * omega / (C_k * temperature)) - 1.)

    
    #is i -> j an absorption?
    elif len(transition_data[( transition_data['Initial Index'] == j) &  (transition_data['Final Index'] == i)]) == 1:
        df = transition_data[( transition_data['Initial Index'] == j) &  (transition_data['Final Index'] == i)]
        omega = 2.0 * np.pi * cts.c/  (df['Wavelength (nm)'].iloc[0]*1e-9)       
        gamma =df['Transition rate (s-1)'].iloc[0]
        
        modeOccupationTerm = 0
        try:
            degeneracyTerm = (2 * df['Initial J'].iloc[0]  +1) / (2*df['Final J'].iloc[0]   + 1)
        except TypeError:
            degeneracyTerm = (2 * float(Fraction(df['Initial J'].iloc[0]))  +1) / (2* float(Fraction(df['Final J'].iloc[0] ))  + 1)
         # only possible by absorbing thermal photons ?
        if (hbar * omega < 100 * C_k * temperature) and (omega > 1e2):
             modeOccupationTerm += 1. / \
                 (np.exp(hbar * omega / (C_k * temperature)) - 1.)
        
 
    
    # i -> j must be too weak to care about, set the transition rate to 0
    else:
        degeneracyTerm = 0
        modeOccupationTerm = 0
        gamma = 0
    return gamma* degeneracyTerm * modeOccupationTerm
   
"""
######################
### the above function was extermely slow due to Pandas boolean search being exteremly slow. Build a pre-indexed dictionary for 50x speedup.
######################
# Build a dict: (Initial, Final) -> row
trans_dict = {
    (row['Initial Index'], row['Final Index']): row
    for _, row in transition_data.iterrows()
} #find all the indeces in advance
def calculate_Transition_Rate(i, j, temperature=0,uncertainty=False):
    key_decay = (i, j)
    key_abs   = (j, i)

    row = trans_dict.get(key_decay)
    
    if uncertainty == False:
    
        if row is not None:
            # decay i -> j
            omega = 2.0 * np.pi * cts.c / (row['Wavelength (nm)'] * 1e-9)
            gamma = row['Transition rate (s-1)']
        
            modeOccupationTerm = 1.0
            degeneracyTerm = 1.0

            if (hbar * omega < 100 * C_k * temperature) and (omega > 1e2):
                modeOccupationTerm += 1.0 / (
                    np.exp(hbar * omega / (C_k * temperature)) - 1.0
                    )

        else:
            # maybe absorption j -> i
            row = trans_dict.get(key_abs)
            if row is not None:
                omega = 2.0 * np.pi * cts.c / (row['Wavelength (nm)'] * 1e-9)
                gamma = row['Transition rate (s-1)']
                
                # ideally preprocess J columns to float once instead of using Fraction
                try:
                    degeneracyTerm = (2 * row['Initial J'] + 1) / (2 * row['Final J'] + 1)
                except TypeError:
                    degeneracyTerm = (
                            2 * float(Fraction(row['Initial J'])) + 1
                            ) / (2 * float(Fraction(row['Final J'])) + 1)

                modeOccupationTerm = 0.0
                if (hbar * omega < 100 * C_k * temperature) and (omega > 1e2):
                    modeOccupationTerm += 1.0 / (
                        np.exp(hbar * omega / (C_k * temperature)) - 1.0
                        )

            else:
                # too weak / non-existent
                gamma = 0.0
                degeneracyTerm = 0.0
                modeOccupationTerm = 0.0
    elif uncertainty == True:
    
        if row is not None:
            # decay i -> j
            omega = 2.0 * np.pi * cts.c / (row['Wavelength (nm)'] * 1e-9)
            gamma = row['Transition rate (s-1)']
            # ADD MONTE CARLO RANDOMNESS
            gamma += np.random.normal(0, row['Transition rate uncertainty'])
            modeOccupationTerm = 1.0
            degeneracyTerm = 1.0

            if (hbar * omega < 100 * C_k * temperature) and (omega > 1e2):
                modeOccupationTerm += 1.0 / (
                    np.exp(hbar * omega / (C_k * temperature)) - 1.0
                    )

        else:
            # maybe absorption j -> i
            row = trans_dict.get(key_abs)
            if row is not None:
                omega = 2.0 * np.pi * cts.c / (row['Wavelength (nm)'] * 1e-9)
                gamma = row['Transition rate (s-1)']
                # ADD MONTE CARLO RANDOMNESS
                gamma += np.random.normal(0, row['Transition rate uncertainty'])
                # ideally preprocess J columns to float once instead of using Fraction
                try:
                    degeneracyTerm = (2 * row['Initial J'] + 1) / (2 * row['Final J'] + 1)
                except TypeError:
                    degeneracyTerm = (
                            2 * float(Fraction(row['Initial J'])) + 1
                            ) / (2 * float(Fraction(row['Final J'])) + 1)

                modeOccupationTerm = 0.0
                if (hbar * omega < 100 * C_k * temperature) and (omega > 1e2):
                    modeOccupationTerm += 1.0 / (
                        np.exp(hbar * omega / (C_k * temperature)) - 1.0
                        )

            else:
                # too weak / non-existent
                gamma = 0.0
                degeneracyTerm = 0.0
                modeOccupationTerm = 0.0

    return gamma * degeneracyTerm * modeOccupationTerm

#%%  
#list the state index for all 2F Term states
#list the state index for all 2D Term states
F_Term_indexes = []
D_Term_indexes = []

for i in range(len(state_data)):
    if state_data.iloc[i]['Term'] == '2F':
        F_Term_indexes.append(i)
        #print(state_data.iloc[i])
    elif state_data.iloc[i]['Term'] == '2D':
        D_Term_indexes.append(i)
        #print(state_data.iloc[i])

print(F_Term_indexes)

print(D_Term_indexes)

#check if there rate is non-zero in arc.
for f in F_Term_indexes:
    for d in D_Term_indexes:
        #print('%r --> %r transiton' %(f,d))
        arc_tuple=(int(state_data.iloc[f]['Configuration'][:-1]),3,float(Fraction(state_data.iloc[f]['J'])), int(state_data.iloc[d]['Configuration'][:-1]),2,float(Fraction(state_data.iloc[d]['J'])))
        try: 
            arc_rate = atom.getTransitionRate(*arc_tuple)
        except ValueError:
            arc_rate = 0
        if arc_rate >= 1:
#if true, check if the rate is incorrectly zero in the rate file
            portal_rate = calculate_Transition_Rate(f, d)
            if portal_rate <=1:
#if the rate is zero in the rate file, append this to the transition list with the  rate given in ARC

                arc_wav = -1*atom.getTransitionWavelength(*arc_tuple)*1e9
                arc_matrix_element = abs(atom.getReducedMatrixElementJ(*arc_tuple))
                print('missing %r --> %r transiton, rate = %r, lambda = %r' %(f,d, arc_rate, arc_wav))
                initial_state = state_data.iloc[f]
                final_state = state_data.iloc[d]

                new_row = pd.DataFrame([{'Initial Configuration': initial_state['Configuration'], 
                           'Initial term': initial_state['Term'],
                           'Initial J': initial_state['J'],
                           'Final configuration': final_state['Configuration'],
                           'Final term': final_state['Term'],
                           'Final J': final_state['J'],
                           'Wavelength (nm)': arc_wav,
                           'Wavelength uncertainty (nm)': np.nan,
                           'Wavelength w/uncertainty (nm)': np.nan,
                           'Matrix element (a.u)': arc_matrix_element,
                           'Matrix el. uncertainty (a.u)': np.nan,
                           'Matrix el. w/uncertainty (a.u)': np.nan,
                           'Transition rate (s-1)': arc_rate,
                           'Transition rate uncertainty': 0.1*arc_rate,
                           'Transition rate w/uncertainty (s-1)':np.nan,
                           'Branching ratio': np.nan,
                           'Branching ratio uncertainty':np.nan,
                           'Branching ratio w/uncertainty': np.nan,
                           'Lifetime reference': 'ARC',
                           'Matrix Element Reference': 'ARC',
                           'Initial Index': f,
                           'Final Index': d}
                           ])
                transition_data = pd.concat([transition_data, new_row],ignore_index = True)

#%%
transition_data.to_excel(new_file,index = False)
