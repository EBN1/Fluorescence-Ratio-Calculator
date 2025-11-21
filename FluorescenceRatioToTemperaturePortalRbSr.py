"""
a script to determine ratios of fluorescence at different wavelengths as a function of temperature
"
@author: ebn
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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
from tabulate import tabulate

plt.style.use('default')
plt.style.use('paper')

species = 'Rb'  #Currently you may choose 'Rb' or 'Sr'

os.chdir(r'C:\Users\ebn1\OneDrive - NIST\Rydberg Thermometry\Fluorescence-Ratio-Calculator')  #Change this to the directory of the .py file

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
   


#%% All the theory stuff and defined functions here



##############
# Define some useful functions
##############
def sortStates(states):  #is this still needed?
    #This method sorts the states being considered by energy levels and returns the adjusted state and energy numpy arrays
    energies = np.zeros((len(states)))
    
    for ii, state in enumerate(states):
        energies[ii] = atom.getEnergy(state['n'], state['L'], state['J'], s=ss)

    # Now sort the states according to energy:
    inds = np.argsort(energies)
    energies = energies[inds]
    states = states[inds]
    
    return states, energies

def fillGammas(states, temp = 298, target_state = 6, Fg = 3 ):
    # A matrix to hold the decay rates and the normalization factor is returned
    GAMMA = np.zeros((len(states), len(states)))
    fails = 0
    #number of bright Zeeman sublevels over total number of Zeeman sublevels.  I can't think of a smart way to get these numbers from the constituent angular momentum, but there are just 4 cases and I can count them
    
    degeneracy_correction = 1
    try:
        if state_data['J'].iloc[target_state] == '3/2' and Fg == 3:
            degeneracy_correction = 18/24
        elif state_data['J'].iloc[target_state] == '1/2' and Fg == 3:
            degeneracy_correction = 11/12
        elif state_data['J'].iloc[target_state] == '3/2' and Fg == 2:
            degeneracy_correction = 12/24
        elif state_data['J'].iloc[target_state] == '1/2' and Fg == 2:
            degeneracy_correction = 9/12
        elif state_data['J'].iloc[target_state] == 1 and Fg == 0:
                degeneracy_correction = 1/3
        elif state_data['J'].iloc[target_state] == 0 and Fg == 1:
                degeneracy_correction = 1
    except:
        fails += 1
        print("%r, %r, %r" % ( state_data['J'].iloc[target_state], Fg, degeneracy_correction))
        print('failed at 1')
    # A nested for loop to build the GAMMA array
    fails = 0
    for ii in range(len(states)):
        for jj in range(len(states)):
            # Make sure the connecting state is higher in energy and separated by
            # one unit of orbital angular momentum (primary selection rule for E1
            # transitions)
            
            
            if jj!=ii :
                #Assign transition rate as normal if spontaneous decay, or stimulated out of something other than the target_state
                if ii!= target_state:
                    try:
                        #gamma = atom.getTransitionRate(y["n"], y["L"], y["J"], x["n"], x["L"], x["J"], s=ss,temperature=temp)
                        gamma = calculate_Transition_Rate(ii,jj, temperature = temp)
                        GAMMA[jj][ii] += gamma
                        GAMMA[ii][ii] -= gamma  #was previously [jj][jj]
                        #print('success at x: %r, %r, %r,   y: %r, %r, %r'  %(x["n"], x["L"], x["J"], y["n"], y["L"], y["J"]))
              
                    
                    except:
                        fails += 1
                        print('failed at 2 x: %r,,   y: '  %(ii, jj))
                else:
                    try:
                        gamma = degeneracy_correction*calculate_Transition_Rate(ii,jj, temperature = temp)
                        GAMMA[jj][ii] += gamma
                        GAMMA[ii][ii] -= gamma  #was previously [jj][jj]
                    except:
                        fails += 1
                        print('failed at 3 x: %r,,   y: '  %(ii, jj))
        
       
        
    #Normalize gamma
    gamma_norm = np.amax(np.abs(GAMMA))
    GAMMA = GAMMA/gamma_norm
    print("Failures:", fails)
    return GAMMA, gamma_norm

def integrate(distribution, GAMMA, gamma_norm, time, max_step=np.inf):
    #Set an initial population array, and then evolve.
    N0 = distribution
    sol = solve_ivp(lambda t, y: GAMMA @ y, [0, np.amax(time)], N0, dense_output = True, max_step=max_step)
        
    print("Done")
    return sol

def findState(parameters, States): #do we still need this?
    state_of_interest = np.array(parameters, dtype=[('n','i4'), ('L', 'i4'), ('J', 'f4')])
    ind = state_of_interest==States
    return ind

def find_equilibrium_population(R, GAMMA):
    # Find the singular values:
    U, S, VH = np.linalg.svd(R+GAMMA)

    Neq = np.compress(S <= 1e-10, VH, axis=0).T
    Neq /= np.sum(Neq)

    if Neq.shape[1] > 1:
        Neq = np.nan*Neq[:, 0]
        #raise ValueError("more than one equilbrium state found")

    # The operations above return a column vector, this collapses the column
    # vector into an array:
    Neq = Neq.flatten()
    
    return Neq

def simulate_spectrum(Neqs, states, transitions, temp = 300):
    wav = np.array(transitions['Wavelength (nm)'])
    intensity = []
    transition_label = np.array(transitions.index) #np.array([transitions['Initial Index'], transitions['Final Index']])
    
    for kk in range(len(transitions)):
        intensity.append(calculate_Transition_Rate(transitions['Initial Index'].iloc[kk], transitions['Final Index'].iloc[kk], temperature=temp)*\
                                                    Neqs[transitions['Initial Index'].iloc[kk]]  )
            
            
    return wav, np.array(intensity), transition_label

def superimpose_gaussians(x, centers, heights, sigma):
    peaks = np.zeros(x.shape)
    
    for center, height in zip(centers, heights):
        peaks += height*np.exp(-(x-center)**2/2/sigma**2)
        
    return peaks

def calculate_fluoresence_ratios( transition_1, transition_2, intensities):
    intensity1= intensities[transition_1]
    intensity2= intensities[transition_2]

    return intensity1/intensity2
#######
#Make a list of states to include in the rate equations
#######
#state_list = [] #Li
#state_list = [(3,2,1.5), (3,2,2.5)] #K
#state_list = [(4,2,1.5), (4,2,2.5), (4,3,2.5), (4,3,3.5)] #Rb, initialized with the 4D & 4F states, the only ones with n =4
#state_list = [(5,2,1.5), (5,2,2.5), (4,3,2.5), (4,3,3.5), (5,3,2.5), (5,3,3.5), (5,4,3.5), (5,4,4.5)] # Cs, initialized with the 4D & 4F states, the only ones with n =4
#state_list = [(4,2,2),(4,3,3)]#Strontium


#%%

########
#Set up the rate equations
########

#laser induced transitions
if species == 'Sr':
    
    targets = [[get_index('5s.6p', '1P', 1),0],
               [get_index('5s.7p', '1P', 1),0],
               [get_index('5s.6s', '1S', 0),1],
               [get_index('5s.6s', '1S', 0),1]]
    
    
    Excitation_Rate = 10  #I believe this is in units of Gamma

    GAMMA, gamma_norm = fillGammas(state_data, temp = 298, target_state= 6, Fg=0)  #this is here so just to get a matrix of the neccessary shape
    
    # Make the pumping matrix
    R = np.zeros(GAMMA.shape)
    #single photon excitation
        
    ind1 = get_index('5s2', '1S', 0) 
    ind2 = targets[0][0]
        

    
    R[ind1, ind1] -= Excitation_Rate
    R[ind1, ind2] += Excitation_Rate
    R[ind2, ind2] -= Excitation_Rate
    R[ind2, ind1] += Excitation_Rate
    
    R1=R  #use this if exciting 5S --> 6P



    # Make the pumping matrix
    R = np.zeros(GAMMA.shape)
    #single photon excitation
    
    ind1 = get_index('5s2', '1S', 0) 
    ind2 = targets[1][0]
    
    R[ind1, ind1] -= Excitation_Rate
    R[ind1, ind2] += Excitation_Rate
    R[ind2, ind2] -= Excitation_Rate
    R[ind2, ind1] += Excitation_Rate

    R2=R  #use this if exciting 5S --> 7P


    # Make the pumping matrix
    R = np.zeros(GAMMA.shape)
    #two photon excitation
    ind1 = get_index('5s2', '1S', 0) 
    ind2 = get_index('5s.5p', '1P', 1)
    ind3 = targets[2][0]
    
    R[ind1, ind1] -= Excitation_Rate
    R[ind1, ind2] += Excitation_Rate
    R[ind2, ind2] -= Excitation_Rate
    R[ind2, ind1] += Excitation_Rate

    R[ind2, ind2] -= Excitation_Rate
    R[ind2, ind3] += Excitation_Rate
    R[ind3, ind3] -= Excitation_Rate
    R[ind3, ind2] += Excitation_Rate

    R3=R  #use this if exciting 5S --> 5P --> 7S



    # Make the pumping matrix
    R = np.zeros(GAMMA.shape)
    #two photon excitation
    ind1 = get_index('5s2', '1S', 0) 
    ind2 = get_index('5s.5p', '1P', 1)
    ind3 = targets[3][0]
    
    R[ind1, ind1] -= Excitation_Rate
    R[ind1, ind2] += Excitation_Rate
    R[ind2, ind2] -= Excitation_Rate
    R[ind2, ind1] += Excitation_Rate

    R[ind2, ind2] -= Excitation_Rate
    R[ind2, ind3] += Excitation_Rate
    R[ind3, ind3] -= Excitation_Rate
    R[ind3, ind2] += Excitation_Rate

    R4=R  #use this if exciting 5S --> 5P --> 8S

    # List of transitions we want to model

    #Rb
    transition_0 = transition_data[( transition_data['Initial Index'] == get_index('5s.6p','1P',1)) &\
                                      (transition_data['Final Index'] == get_index('5s2', '1S', 0))].index[0]
    transition_1 = transition_data[( transition_data['Initial Index'] == get_index('5s.5p','1P',1)) &\
                                      (transition_data['Final Index'] == get_index('5s2', '1S', 0))].index[0]
    transition_2 = transition_data[( transition_data['Initial Index'] == get_index('5s.7p','1P',1)) &\
                                      (transition_data['Final Index'] == get_index('5s2', '1S', 0))].index[0]


    transition_3 = transition_data[( transition_data['Initial Index'] == get_index('5s.5d','1D',2)) &\
                                      (transition_data['Final Index'] == get_index('5s.5p', '1P', 1))].index[0]

    transition_4 = transition_data[( transition_data['Initial Index'] == get_index('5s.6d','1D',2)) &\
                                      (transition_data['Final Index'] == get_index('5s.5p', '1P', 1))].index[0]


    transition_5 = transition_data[( transition_data['Initial Index'] == get_index('5s.6s','1S',0)) &\
                                      (transition_data['Final Index'] == get_index('5s.5p', '1P', 1))].index[0]
    
    transition_6 = transition_data[( transition_data['Initial Index'] == get_index('5s.7s','1S',0)) &\
                                      (transition_data['Final Index'] == get_index('5s.5p', '1P', 1))].index[0]

#    transition_7 = transition_data[( transition_data['Initial Index'] == get_index('8','1S',0)) &\
#                                      (transition_data['Final Index'] == get_index('5p', '1P', 1))].index[0]
#    transition_8 = transition_data[( transition_data['Initial Index'] == get_index('8','1S',0)) &\
#                                      (transition_data['Final Index'] == get_index('6p', '1P', 1))].index[0]
#    transition_9 = transition_data[( transition_data['Initial Index'] == get_index('8','1S',0)) &\
#                                      (transition_data['Final Index'] == get_index('7p', '1P', 1))].index[0]
#    transition_10 = transition_data[( transition_data['Initial Index'] == get_index('9','1S',0)) &\
#                                      (transition_data['Final Index'] == get_index('5p', '1P', 1))].index[0]
 


#laser induced transitions
if species == 'Rb':
    
    targets = [[get_index('6p', '2P', '1/2'),3],
               [get_index('6p', '2P', '3/2'),3],
               [get_index('7p', '2P', '1/2'),3],
               [get_index('7p', '2P', '3/2'),3]]
    
    Excitation_Rate = 10  #I believe this is in units of Gamma

    GAMMA, gamma_norm = fillGammas(state_data, temp = 298, target_state= 6, Fg=0) #this is here so just to get a matrix of the neccessary shape
    
    # Make the pumping matrix
    R = np.zeros(GAMMA.shape)
    #single photon excitation

    
    ind1 = get_index('5s', '2S','1/2') 
    ind2 = targets[0][0]
    
    
    R[ind1, ind1] -= Excitation_Rate
    R[ind1, ind2] += Excitation_Rate
    R[ind2, ind2] -= Excitation_Rate
    R[ind2, ind1] += Excitation_Rate
    
    R1=R  #use this if exciting 5S --> 6P



    # Make the pumping matrix
    R = np.zeros(GAMMA.shape)
    #single photon excitation
    ind1 = get_index('5s', '2S','1/2') 
    ind2 = targets[1][0]
    
    R[ind1, ind1] -= Excitation_Rate
    R[ind1, ind2] += Excitation_Rate
    R[ind2, ind2] -= Excitation_Rate
    R[ind2, ind1] += Excitation_Rate

    R2=R  #use this if exciting 5S --> 7P


    # Make the pumping matrix
    R = np.zeros(GAMMA.shape)
    #two photon excitation
    ind1 = get_index('5s', '2S','1/2') 
    ind2 = targets[2][0]
    
    R[ind1, ind1] -= Excitation_Rate
    R[ind1, ind2] += Excitation_Rate
    R[ind2, ind2] -= Excitation_Rate
    R[ind2, ind1] += Excitation_Rate


    R3=R  #use this if exciting 5S --> 5P --> 7S



    # Make the pumping matrix
    R = np.zeros(GAMMA.shape)
    ind1 = get_index('5s', '2S','1/2') 
    ind2 = targets[3][0]
    
    R[ind1, ind1] -= Excitation_Rate
    R[ind1, ind2] += Excitation_Rate
    R[ind2, ind2] -= Excitation_Rate
    R[ind2, ind1] += Excitation_Rate

 

    R4=R  #use this if exciting 5S --> 5P --> 8S

    # List of transitions we want to model

    #Rb
    transition_0 =transition_data[( transition_data['Initial Index'] == get_index('5p','2P','1/2')) &\
                                  (transition_data['Final Index'] == get_index('5s', '2S', '1/2'))].index[0]
    transition_1 = transition_data[( transition_data['Initial Index'] == get_index('5p','2P','3/2')) &\
                                  (transition_data['Final Index'] == get_index('5s', '2S', '1/2'))].index[0]
         
        
    transition_2 =transition_data[( transition_data['Initial Index'] == get_index('5d','2D','5/2')) &\
                                  (transition_data['Final Index'] == get_index('5p', '2P', '3/2'))].index[0]
        
    transition_3 =transition_data[( transition_data['Initial Index'] == get_index('5d','2D','3/2')) &\
                                  (transition_data['Final Index'] == get_index('5p', '2P', '3/2'))].index[0]

    transition_4 =transition_data[( transition_data['Initial Index'] == get_index('5d','2D','3/2')) &\
                                  (transition_data['Final Index'] == get_index('5p', '2P', '1/2'))].index[0]

    transition_5 =transition_data[( transition_data['Initial Index'] == get_index('7s','2S','1/2')) &\
                                  (transition_data['Final Index'] == get_index('5p', '2P', '3/2'))].index[0]
    transition_6 =transition_data[( transition_data['Initial Index'] == get_index('7s','2S','1/2')) &\
                                  (transition_data['Final Index'] == get_index('5p', '2P', '1/2'))].index[0]    
        

    transition_7 =transition_data[( transition_data['Initial Index'] == get_index('6d','2D','5/2')) &\
                                  (transition_data['Final Index'] == get_index('5p', '2P', '3/2'))].index[0]
        
    transition_8 =transition_data[( transition_data['Initial Index'] == get_index('6d','2D','3/2')) &\
                                  (transition_data['Final Index'] == get_index('5p', '2P', '3/2'))].index[0]

    transition_9 =transition_data[( transition_data['Initial Index'] == get_index('6d','2D','3/2')) &\
                                  (transition_data['Final Index'] == get_index('5p', '2P', '1/2'))].index[0]
    
    transition_10 =transition_data[( transition_data['Initial Index'] == get_index('8s','2S','1/2')) &\
                                  (transition_data['Final Index'] == get_index('5p', '2P', '3/2'))].index[0]
    transition_11 =transition_data[( transition_data['Initial Index'] == get_index('8s','2S','1/2')) &\
                                  (transition_data['Final Index'] == get_index('5p', '2P', '1/2'))].index[0]  
    
    
    transition_12 =transition_data[( transition_data['Initial Index'] == get_index('7p','2P','1/2')) &\
                                  (transition_data['Final Index'] == get_index('5s', '2S', '1/2'))].index[0]
    transition_13 = transition_data[( transition_data['Initial Index'] == get_index('7p','2P','3/2')) &\
                                  (transition_data['Final Index'] == get_index('5s', '2S', '1/2'))].index[0]
    
    transition_14 =transition_data[( transition_data['Initial Index'] == get_index('8p','2P','1/2')) &\
                                 (transition_data['Final Index'] == get_index('5s', '2S', '1/2'))].index[0]
    transition_15 = transition_data[( transition_data['Initial Index'] == get_index('8p','2P','3/2')) &\
                                 (transition_data['Final Index'] == get_index('5s', '2S', '1/2'))].index[0]



#Ts = np.sort(np.unique(Approximate_Temperatures))
#Ts = np.sort(np.unique(np.concatenate([Approximate_Temperatures,[273,373] ])))
#Ts_calibration_index = np.where( Ts == T_calibrate)[0][0]

Ts= np.arange(230,381,25)
ratios = []
ratios1 = []
ratios2 = []
ratios3 = []
ratios4=[]


if species == 'Rb':
    for T in Ts:
        # Remake gamma matrix at this temperature:
    

    # Make the decay matrix:
        GAMMA, gamma_norm = fillGammas(state_data, temp = T, target_state= targets[0][0], Fg=targets[0][1])
        Neq1 = find_equilibrium_population(R1, GAMMA)
        
    # Make the decay matrix:
        GAMMA, gamma_norm = fillGammas(state_data, temp = T, target_state= targets[1][0], Fg=targets[1][1])
        Neq2 = find_equilibrium_population(R2, GAMMA)
    
    # Make the decay matrix:
        GAMMA, gamma_norm = fillGammas(state_data, temp = T, target_state= targets[2][0], Fg=targets[2][1])
        Neq3 = find_equilibrium_population(R3, GAMMA)
    
    # Make the decay matrix:
        GAMMA, gamma_norm = fillGammas(state_data, temp = T, target_state= targets[3][0], Fg=targets[3][1])
        Neq4 = find_equilibrium_population(R4, GAMMA)

    # Simulate the spectrum:

        wavs1, intensities1, transitions1 = simulate_spectrum(Neq1, state_data, transition_data, temp = T)
        wavs2, intensities2, transitions2 = simulate_spectrum(Neq2, state_data, transition_data, temp = T)
        wavs3, intensities3, transitions3 = simulate_spectrum(Neq3, state_data, transition_data, temp = T)
        wavs4, intensities4, transitions4 = simulate_spectrum(Neq4, state_data, transition_data, temp = T)

    # Find the fluorescence ratios:
        ratios1.append([0,
                   calculate_fluoresence_ratios(transition_0, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_2, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_3, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_4, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_5, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_6, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_7, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_8, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_9, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_10, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_11, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_12, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_13, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_14, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_15, transition_1, intensities1),
                   
                  ])
    
        ratios2.append([0,
                   calculate_fluoresence_ratios(transition_0, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_2, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_3, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_4, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_5, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_6, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_7, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_8, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_9, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_10, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_11, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_12, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_13, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_14, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_15, transition_1, intensities2),
                   
                  ])
        ratios3.append([0,
                   calculate_fluoresence_ratios(transition_0, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_2, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_3, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_4, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_5, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_6, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_7, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_8, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_9, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_10, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_11, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_12, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_13, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_14, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_15, transition_1, intensities3),
                   
                  ])
        ratios4.append([0,
                   calculate_fluoresence_ratios(transition_0, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_2, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_3, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_4, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_5, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_6, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_7, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_8, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_9, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_10, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_11, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_12, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_13, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_14, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_15, transition_1, intensities4),
                   
                  ])
elif species == 'Sr':
    for T in Ts:
        # Remake gamma matrix at this temperature:
    

    # Make the decay matrix:
        GAMMA, gamma_norm = fillGammas(state_data, temp = T, target_state= targets[0][0], Fg=targets[0][1])
        Neq1 = find_equilibrium_population(R1, GAMMA)
        
    # Make the decay matrix:
        GAMMA, gamma_norm = fillGammas(state_data, temp = T, target_state= targets[1][0], Fg=targets[1][1])
        Neq2 = find_equilibrium_population(R2, GAMMA)
    
    # Make the decay matrix:
        GAMMA, gamma_norm = fillGammas(state_data, temp = T, target_state= targets[2][0], Fg=targets[2][1])
        Neq3 = find_equilibrium_population(R3, GAMMA)
    
    # Make the decay matrix:
        GAMMA, gamma_norm = fillGammas(state_data, temp = T, target_state= targets[3][0], Fg=targets[3][1])
        Neq4 = find_equilibrium_population(R4, GAMMA)

    # Simulate the spectrum:

        wavs1, intensities1, transitions1 = simulate_spectrum(Neq1, state_data, transition_data, temp = T)
        wavs2, intensities2, transitions2 = simulate_spectrum(Neq2, state_data, transition_data, temp = T)
        wavs3, intensities3, transitions3 = simulate_spectrum(Neq3, state_data, transition_data, temp = T)
        wavs4, intensities4, transitions4 = simulate_spectrum(Neq4, state_data, transition_data, temp = T)

    # Find the fluorescence ratios:
        ratios1.append([0,
                   calculate_fluoresence_ratios(transition_0, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_2, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_3, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_4, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_5, transition_1, intensities1),
                   calculate_fluoresence_ratios(transition_6, transition_1, intensities1),
                   
                   
                  ])
    
        ratios2.append([0,
                   calculate_fluoresence_ratios(transition_0, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_2, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_3, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_4, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_5, transition_1, intensities2),
                   calculate_fluoresence_ratios(transition_6, transition_1, intensities2),
                   
                   
                  ])
        ratios3.append([0,
                   calculate_fluoresence_ratios(transition_0, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_2, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_3, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_4, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_5, transition_1, intensities3),
                   calculate_fluoresence_ratios(transition_6, transition_1, intensities3),
                  
                  ])
        ratios4.append([0,
                   calculate_fluoresence_ratios(transition_0, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_2, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_3, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_4, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_5, transition_1, intensities4),
                   calculate_fluoresence_ratios(transition_6, transition_1, intensities4),
                   
                  ])
 
ratios1 = np.array(ratios1)
ratios2 = np.array(ratios2)
ratios3 = np.array(ratios3)
ratios4 = np.array(ratios4)

#%%
#############
# Print some stuff to prove we caluclated some ratios
#############
Ts_calibration_index=2  #at which temperature in the array Ts would you like to compare ratios?
if species =='Rb':
    transition = np.array([transition_0, transition_1, transition_2, transition_3, transition_4, transition_5, transition_6, transition_7, transition_8, transition_9, transition_10, transition_11, transition_12, transition_13, transition_14, transition_15])
    header = ["i",
                                   "$\lambda$ (nm)",
                                   "Config_i", 
                                   "Config f",
                                   "6P_1/2 Signal @ T= %d K"%Ts[Ts_calibration_index], "6P_3/2 Signal @ T=%d K"%Ts[Ts_calibration_index],
                                   "7P_1/2 Signal @ T=%d K"%Ts[Ts_calibration_index], "7P_3/2 Signal @ T=%d K"%Ts[Ts_calibration_index] ]
elif species == 'Sr':
    transition = np.array([transition_0, transition_1, transition_2, transition_3, transition_4, transition_5, transition_6])
    header = ["i",
                                   "$\lambda$ (nm)",
                                   "Config_i", 
                                   "Config f",
                                   "6P_1 Signal @ T= %d K"%Ts[Ts_calibration_index], "7P_1 Signal @ T=%d K"%Ts[Ts_calibration_index],
                                   "7S_0 Signal @ T=%d K"%Ts[Ts_calibration_index], "8S_0 Signal @ T=%d K"%Ts[Ts_calibration_index] ]

table = []
for j, i in enumerate(transition):
    if i == 0:
        table.append([i, transition_data['Wavelength (nm)'].iloc[i],
                      transition_data['Initial Configuration'].iloc[i]+
                      transition_data['Initial term'].iloc[i]+
                      str(transition_data['Initial J'].iloc[i]),
                      transition_data['Final configuration'].iloc[i]+ 
                      transition_data['Final term'].iloc[i]+ 
                      str(transition_data['Final J'].iloc[i]),
                 ratios1[Ts_calibration_index,1], ratios2[Ts_calibration_index,1], ratios3[Ts_calibration_index,1], ratios4[Ts_calibration_index,1]])
    
    elif i == 1:
        table.append([i,  transition_data['Wavelength (nm)'].iloc[i],
                     transition_data['Initial Configuration'].iloc[i]+
                     transition_data['Initial term'].iloc[i]+
                     str(transition_data['Initial J'].iloc[i]),
                     transition_data['Final configuration'].iloc[i]+ 
                     transition_data['Final term'].iloc[i]+ 
                     str(transition_data['Final J'].iloc[i]),
                      1, 1,1,1])
    else:
        table.append([i,  transition_data['Wavelength (nm)'].iloc[i],
                      transition_data['Initial Configuration'].iloc[i]+
                      transition_data['Initial term'].iloc[i]+
                      str(transition_data['Initial J'].iloc[i]),
                      transition_data['Final configuration'].iloc[i]+ 
                      transition_data['Final term'].iloc[i]+ 
                      str(transition_data['Final J'].iloc[i]),
                      ratios1[Ts_calibration_index,j], ratios2[Ts_calibration_index,j],ratios3[Ts_calibration_index,j], ratios4[Ts_calibration_index,j]])
    
    #print("%d, %r, %r" %  (i, transition[i], -atom.getTransitionWavelength(int(transition[i,0]),int(transition[i,1]),transition[i,2],int(transition[i,3]),int(transition[i,4]),transition[i,5])
    #                      ))
    


#%%    


print(tabulate(table,headers= header))


#####
# Calculate the expected ratio for the filters

#Filter number   #Filter Wavelength
#1                  780
#2                  630
#3                  760
#4                  740
#5                  850
#6                  620
#######

if species == 'Sr':

    Predicted_ratios_2_1 = np.array([ratios1[:,1], ratios2[:,1], ratios3[:,1], ratios4[:,1]])

    Predicted_ratios_3_1 = np.array([ratios1[:,2], ratios2[:,2], ratios3[:,2], ratios4[:,2]])

    Predicted_ratios_4_1 = np.array([ratios1[:,3], ratios2[:,3], ratios3[:,3], ratios4[:,3]])

    Predicted_ratios_5_1 = np.array([ratios1[:,4], ratios2[:,4], ratios3[:,4], ratios4[:,4]])

    Predicted_ratios_6_1 = np.array([ratios1[:,5], ratios2[:,5], ratios3[:,5], ratios4[:,5]])
    
if species == 'Rb':

    Predicted_ratios_2_1 = np.array([ratios1[:,7]+ratios1[:,8], ratios2[:,7]+ratios2[:,8],ratios3[:,7]+ratios3[:,8], ratios4[:,7]+ratios4[:,8] ])

    Predicted_ratios_3_1 = np.array([ratios1[:,4], ratios2[:,4], ratios3[:,4], ratios4[:,4] ])

    Predicted_ratios_4_1 = np.array([ratios1[:,5], ratios2[:,5], ratios3[:,5], ratios4[:,5] ])

    Predicted_ratios_5_1 = np.array([ratios1[:,11], ratios2[:,11], ratios3[:,11], ratios4[:,11]])

    Predicted_ratios_6_1 = np.array([ratios1[:,9]+ratios1[:,10], ratios2[:,9]+ratios2[:,10],ratios3[:,9]+ratios3[:,10], ratios4[:,9]+ratios4[:,10] ])


Predicted_ratios_2_3 = Predicted_ratios_2_1/Predicted_ratios_3_1

Predicted_ratios_2_4 =Predicted_ratios_2_1/Predicted_ratios_4_1

Predicted_ratios_3_4 =Predicted_ratios_3_1/Predicted_ratios_4_1

Predicted_ratios_4_3 = Predicted_ratios_2_3/Predicted_ratios_2_4

#Predicted_ratios_2_1 = np.array(table[7][-4:])+np.array(table[8][-4:])
#Predicted_ratios_6_1 = np.array(table[9][-4:])+np.array(table[10][-4:])

#Predicted_ratios_2_6 = Calibration_ratios_2_1/Calibration_ratios_6_1




#np.save('Predicted_ratios_2_3_portal_degeneracy_correction_ratio_down.npy',Predicted_ratios_2_3)
#np.save('Predicted_ratios_2_4_portal_degeneracy_correction_ratio_down.npy',Predicted_ratios_2_4)
#np.save('Predicted_ratios_3_4_portal_degeneracy_correction_ratio_down.npy',Predicted_ratios_3_4)


#%%
#%%
plt.figure()
plt.plot(Ts, Predicted_ratios_4_3[1])
plt.xlabel("T (K)")
plt.ylabel('ratio 740 nm/760 nm')
plt.title('Rb 6P excitation')


plt.figure()
plt.plot(Ts, Predicted_ratios_4_3[3])
plt.xlabel("T (K)")
plt.ylabel('ratio 740 nm/760 nm')
plt.ylim(0,15)
plt.title('Rb 7P excitation')
#%%

gs_kw = dict(height_ratios=[1.5,1], wspace=.25, hspace=0.3, top =.9, left = .12)
fig, ax = plt.subplots(2,3,figsize=(6.75,3.25) ,  gridspec_kw= gs_kw)
fig.markersize=12
#fig.set_figheight(3.25*2)
#fig.set_figwidth(6.75)



#%% A potentially useful plot for Rb


ax[0,0].plot(Ts, Predicted_ratios_2_1[0],color='pink' , label="$6^1P_1$")
ax[0,0].semilogy(Ts, Predicted_ratios_2_1[1],color='b')
ax[0,0].plot(Ts, Predicted_ratios_2_1[2],color='orange')
ax[0,0].semilogy(Ts, Predicted_ratios_2_1[3],color='green')
ax[0,0].set_title("%d nm/%d nm" %(table[7][1], table[1][1]))

ax[0,1].plot(Ts,    Predicted_ratios_3_1[0],color='pink')
ax[0,1].semilogy(Ts, Predicted_ratios_3_1[1],color='b')
ax[0,1].plot(Ts,    Predicted_ratios_3_1[2],color='orange')
ax[0,1].semilogy(Ts, Predicted_ratios_3_1[3],color='green')
ax[0,1].set_title("%d nm/%d nm" %(table[4][1], table[1][1]))

ax[0,2].plot(Ts,    Predicted_ratios_4_1[0],color='pink')
ax[0,2].semilogy(Ts, Predicted_ratios_4_1[1],color='b')
ax[0,2].plot(Ts,    Predicted_ratios_4_1[2],color='orange')
ax[0,2].semilogy(Ts, Predicted_ratios_4_1[3],color='green')
ax[0,2].set_title("%d nm/%d nm" %(table[5][1], table[1][1]))

ax[1,0].plot(Ts,    Predicted_ratios_5_1[0],color='pink')
ax[1,0].semilogy(Ts, Predicted_ratios_5_1[1],color='b')
ax[1,0].plot(Ts,    Predicted_ratios_5_1[2],color='orange')
ax[1,0].semilogy(Ts, Predicted_ratios_5_1[3],color='green')
ax[1,0].set_title("%d nm/%d nm" %(table[11][1], table[1][1]))


ax[1,1].plot(Ts, Predicted_ratios_6_1[0],color='pink')
ax[1,1].semilogy(Ts, Predicted_ratios_6_1[1],color='b')
ax[1,1].plot(Ts, Predicted_ratios_6_1[2],color='orange')
ax[1,1].semilogy(Ts, Predicted_ratios_6_1[3],color='green')
ax[1,1].set_title("%d nm/%d nm" %(table[9][1], table[1][1]))

ax[1,0].set_xlabel('Temperature (K)')
ax[1,1].set_xlabel('Temperature (K)')
ax[1,2].set_xlabel('Temperature (K)')


ax[1,2].plot(Ts, Predicted_ratios_2_1[0],color='pink' , label="$6^2P_{1/2}$")
ax[1,2].plot(Ts, Predicted_ratios_2_1[1],color='b', label="$6^2P_{3/2}$")
ax[1,2].plot(Ts, Predicted_ratios_2_1[2],color='orange', label="$7^2P_{1/2}$")
ax[1,2].plot(Ts, Predicted_ratios_2_1[3],color='green', label="$7^2P_{3/2}$")
ax[1,2].set_ylim(1,2)
ax[1,2].legend()


#%% Rb summary

gs_kw = dict(width_ratios=[1,1,1,.75], wspace=.25, hspace=0.3, top =.9, left = .12)
fig, ax = plt.subplots(1,4,figsize=(6.75,3.25) , gridspec_kw= gs_kw)
fig.markersize=12
#fig.set_figheight(3.25*2)
#fig.set_figwidth(6.75)



ax[1].plot(Ts, Predicted_ratios_2_3[0],color='pink' , label="$6^1P_1$")
ax[1].semilogy(Ts, Predicted_ratios_2_3[1],color='b')
ax[1].plot(Ts, Predicted_ratios_2_3[2],color='orange')
ax[1].semilogy(Ts, Predicted_ratios_2_3[3],color='green')
ax[1].set_title("%d nm/%d nm" %(table[7][1], table[4][1]))

ax[0].plot(Ts,    Predicted_ratios_2_4[0],color='pink')
ax[0].semilogy(Ts, Predicted_ratios_2_4[1],color='b')
ax[0].plot(Ts,    Predicted_ratios_2_4[2],color='orange')
ax[0].semilogy(Ts, Predicted_ratios_2_4[3],color='green')
ax[0].set_title("%d nm/%d nm" %(table[7][1], table[5][1]))

ax[2].plot(Ts,    Predicted_ratios_3_4[0],color='pink')
ax[2].semilogy(Ts, Predicted_ratios_3_4[1],color='b')
ax[2].plot(Ts,    Predicted_ratios_3_4[2],color='orange')
ax[2].semilogy(Ts, Predicted_ratios_3_4[3],color='green')
ax[2].set_title("%d nm/%d nm" %(table[4][1], table[5][1]))

ax[3].plot(Ts, Predicted_ratios_2_1[0],color='pink' , label="$6^2P_{1/2}$")
ax[3].plot(Ts, Predicted_ratios_2_1[1],color='b', label="$6^2P_{3/2}$")
ax[3].plot(Ts, Predicted_ratios_2_1[2],color='orange', label="$7^2P_{1/2}$")
ax[3].plot(Ts, Predicted_ratios_2_1[3],color='green', label="$7^2P_{3/2}$")
ax[3].set_ylim(1,2)
ax[3].legend()

#%% A potentially useful plot for Sr


ax[0,0].plot(Ts, Predicted_ratios_2_1[0],color='pink' , label="$6^1P_1$")
ax[0,0].semilogy(Ts, Predicted_ratios_2_1[1],color='b')
ax[0,0].plot(Ts, Predicted_ratios_2_1[2],color='orange')
ax[0,0].semilogy(Ts, Predicted_ratios_2_1[3],color='green')
ax[0,0].set_title("%d nm/%d nm" %(table[0][1], table[1][1]))

ax[0,1].plot(Ts,    Predicted_ratios_3_1[0],color='pink')
ax[0,1].semilogy(Ts, Predicted_ratios_3_1[1],color='b')
ax[0,1].plot(Ts,    Predicted_ratios_3_1[2],color='orange')
ax[0,1].semilogy(Ts, Predicted_ratios_3_1[3],color='green')
ax[0,1].set_title("%d nm/%d nm" %(table[2][1], table[1][1]))

ax[0,2].plot(Ts,    Predicted_ratios_4_1[0],color='pink')
ax[0,2].semilogy(Ts, Predicted_ratios_4_1[1],color='b')
ax[0,2].plot(Ts,    Predicted_ratios_4_1[2],color='orange')
ax[0,2].semilogy(Ts, Predicted_ratios_4_1[3],color='green')
ax[0,2].set_title("%d nm/%d nm" %(table[3][1], table[1][1]))

ax[1,0].plot(Ts,    Predicted_ratios_5_1[0],color='pink')
ax[1,0].semilogy(Ts, Predicted_ratios_5_1[1],color='b')
ax[1,0].plot(Ts,    Predicted_ratios_5_1[2],color='orange')
ax[1,0].semilogy(Ts, Predicted_ratios_5_1[3],color='green')
ax[1,0].set_title("%d nm/%d nm" %(table[4][1], table[1][1]))


ax[1,1].plot(Ts, Predicted_ratios_6_1[0],color='pink')
ax[1,1].semilogy(Ts, Predicted_ratios_6_1[1],color='b')

ax[1,1].plot(Ts, Predicted_ratios_6_1[2],color='orange')
ax[1,1].semilogy(Ts, Predicted_ratios_6_1[3],color='green')
ax[1,1].set_title("%d nm/%d nm" %(table[5][1], table[1][1]))

ax[1,0].set_xlabel('Temperature (K)')
ax[1,1].set_xlabel('Temperature (K)')
ax[1,2].set_xlabel('Temperature (K)')


ax[1,2].plot(Ts, Predicted_ratios_2_1[0],color='pink' , label="$6^1P_1$")
ax[1,2].plot(Ts, Predicted_ratios_2_1[1],color='b', label="$7^1P_1$")
ax[1,2].plot(Ts, Predicted_ratios_2_1[2],color='orange', label="$7^1S_0$")
ax[1,2].plot(Ts, Predicted_ratios_2_1[3],color='green', label="$8^1S_0$")
ax[1,2].set_ylim(1,2)
ax[1,2].legend()

#%%

#ax[0,0].fill_between(Ts,Predicted_ratios_2_4_up[3],Predicted_ratios_2_4_down[3], edgecolor='0.5', facecolor ='0.9')

ax[0,0].plot(Ts,Predicted_ratios_2_4[3], color = 'k')
ax[0,0].errorbar(Approximate_Temperatures[number_of_cals:],unumpy.nominal_values(data_2_4[number_of_cals:]), yerr= unumpy.std_devs(data_2_4[number_of_cals:]), color = 'C0', ls='none', marker = '.',markersize =ms)
ax[1,0].set_xlabel('Temperature (K)')
ax[0,0].set_title(r'$r_{\rm{630 nm, 740 nm}}^{(7\,^2\rm{P}_{3/2})}$')
ax[0,0].set_ylabel('Signal Ratio')


ax[0,0].set_xlim(280,350)
ax[0,1].set_xlim(280,350)
ax[0,2].set_xlim(280,350)
ax[1,0].set_xlim(280,350)
ax[1,1].set_xlim(280,350)
ax[1,2].set_xlim(280,350)

ax[0,0].set_ylim(.01,.026)
ax[0,1].set_ylim(.06,0.25)
ax[0,2].set_ylim(0.103,0.163)

#plt.figure('Ratio of 630nm/760 nm, geometric mean of PMT signals, P3/2  to P1/2 calibration')

ax[0,1].fill_between(Ts,Predicted_ratios_2_3_up[3],Predicted_ratios_2_3_down[3], edgecolor='0.5', facecolor ='0.9')

ax[0,1].fill_between(Ts, Predicted_ratios_2_3_up[3]*Predicted_ratios_2_3_up[2]/Predicted_ratios_2_3[2], Predicted_ratios_2_3_down[3]*Predicted_ratios_2_3_down[2]/Predicted_ratios_2_3[2], edgecolor='0.5', facecolor ='0.9')

ax[0,1].plot(Ts,Predicted_ratios_2_3[3], color = 'k')
#ax[0,1].plot(Ts,Predicted_ratios_2_3[2], '-.', color = 'k')
ax[0,1].errorbar(Approximate_Temperatures[number_of_cals:],unumpy.nominal_values(data_2_3[number_of_cals:]), yerr= unumpy.std_devs(data_2_3[number_of_cals:]), color = 'C0', ls='none', marker = '.',markersize =ms)
ax[1,1].set_xlabel('Temperature (K)')
ax[0,1].set_title(r'$r_{\rm{630 nm, 760 nm}}^{(7\,^2\rm{P}_{3/2})}$')


#plt.figure('Ratio of 760 nm/740, geometric mean of PMT signals, P3/2  to P1/2 calibration')

ax[0,2].fill_between(Ts,Predicted_ratios_3_4_up[3],Predicted_ratios_3_4_down[3], edgecolor='0.5', facecolor ='0.9')
ax[0,2].plot(Ts,Predicted_ratios_3_4[3], color = 'k')

ax[0,2].errorbar(Approximate_Temperatures[number_of_cals:],unumpy.nominal_values(data_3_4[number_of_cals:]), yerr= unumpy.std_devs(data_3_4[number_of_cals:]), color = 'C0', ls='none', marker = '.',markersize =ms)
ax[1,2].set_xlabel('Temperature (K)')
ax[0,2].set_title(r'$r_{\rm{760 nm, 740 nm}}^{(7\,^2\rm{P}_{3/2})}$')



