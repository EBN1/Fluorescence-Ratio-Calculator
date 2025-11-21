# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 07:18:30 2025

@author: ebn1

Input:  a transition rate file from the UD atomic portal.  On the portal, select an element>> Transition Rates >> All rates >> Download the xlsx file of transition rates.
Output a list of states extracted from the transition rate file
"""
import numpy as np
import pandas as pd
import os

save_output_file = True
species = 'Rb'


os.chdir(r'C:\Users\cavslab\OneDrive - NIST\CoBRAS\fluorescence\UD Portal files for ARC')
#edit the location of the xlsx file containing the transition rates for the desired species
transition_file ={
    'Sr':  'Sr1-transition-rates.xlsx',
    'Rb':  'Rb1-transition-rates.xlsx',
    }


headers = ['Initial Configuration',
           'Initial term',
           'Initial J',
           'Final configuration',
           'Final term',
           'Final J',
           'Wavelength (nm)',
           'Wavelength error (nm)',
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
           'Wavelength error (nm)': 'float',
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

print('successfully loaded the transition rate file')


#input transition data
#output data frame of states:  index, config, term, J
def get_States(td):
    #states =pd.DataFrame(columns = ['Configuration', 'Term', 'J'])
    df = td[['Initial Configuration',
               'Initial term',
               'Initial J',
               'Final configuration',
               'Final term',
               'Final J']]
    
    states = pd.DataFrame({'Configuration': [df['Initial Configuration'].iloc[0],df['Final configuration'].iloc[0]],
                       'Term': [df['Initial term'].iloc[0],df['Final term'].iloc[0] ],
                       'J': [ df['Initial J'].iloc[0], df['Final J'].iloc[0]] })

    
    
    for i in range(1,len(df)):
        print('states:')
        print(states)
        initial_row = pd.Series({'Configuration': df['Initial Configuration'].iloc[i],
                       'Term': df['Initial term'].iloc[i],
                       'J': df['Initial J'].iloc[i]})
        initial_df = pd.DataFrame({
            'Configuration': [df['Initial Configuration'].iloc[i] ],
            'Term': [df['Initial term'].iloc[i]],
            'J': [df['Initial J'].iloc[i]]
            })
        
        #print('initial row:')
        #print(initial_row)
        #print(initial_df)
        final_row= pd.Series({'Configuration': df['Final configuration'].iloc[i],
                       'Term': df['Final term'].iloc[i],
                       'J': df['Final J'].iloc[i]})
        final_df= pd.DataFrame({
            'Configuration': [df['Final configuration'].iloc[i]],
            'Term': [df['Final term'].iloc[i] ],
            'J': [df['Final J'].iloc[i]]
            })
        
        #print('final row:')
        #print(final_row)
        #print(final_df)
        row_exists = states.apply(lambda row: row.equals(initial_row), axis=1).any()

        #print(f"Does the initial row exist in the DataFrame? {row_exists}")
        
        
        # Check if any row matches the mask
        if not row_exists:
            states= pd.concat([states, initial_df], ignore_index=True)
            
        row_exists = states.apply(lambda row: row.equals(final_row), axis=1).any()
        #print(f"Does the final row exist in the DataFrame? {row_exists}")
        # Check if any row matches the mask
        if not row_exists:
            states= pd.concat( [ states, final_df ], ignore_index=True)
            
    return states 
        # Concatenate the new row to the existing DataFrame
#        df = pd.concat([df, new_row_df], ignore_index=True)


states_data = get_States(transition_data)        

#%%
if save_output_file:
    states_data.to_excel(species + '-states.xlsx')

#%%
import arc
from arc import *
from scipy import constants as cts
from fractions import Fraction


def get_Transition_Rate(i):  #a function to check consistency.  Takes in transition_data index idx, calculated the transition rate from the matrix element and wavelength
    wav = transition_data['Wavelength (nm)'].iloc[i] *1e-9
    print(wav)
    mu  = transition_data['Matrix element (a.u)'].iloc[i]*(cts.e*cts.value('Bohr radius'))
    print(mu)
    Ji = transition_data['Initial J'].iloc[i]
    if isinstance(Ji, str):
        Ji = float(Fraction(Ji))
    print(Ji)
    Jf = transition_data['Final J'].iloc[i]
    if isinstance(Jf, str):
        Jf = float(Fraction(Jf))
    print(Jf)
    Rate_portal = transition_data['Transition rate (s-1)'].iloc[i]
    print(Rate_portal)
    
    rate = (16*np.pi**3 * mu**2)/(3 * cts.epsilon_0*cts.h * wav**3)/(2*Ji+1)
    print(rate)
    return rate, Rate_portal, rate/Rate_portal, Ji, Jf

#atom = arc.Strontium88()
#arc_rate = atom.getTransitionRate(5,1,1,5,0,0,temperature = 0, s=0)

atom = arc.Rubidium()
arc_rate = atom.getTransitionRate(5,1,1.5,5,0,0.5)