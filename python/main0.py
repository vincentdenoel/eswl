import os
import numpy as np
import scipy.io as sio
from eswl import do_buffeting_analysis

# Flag to recalculate the analysis
redo_analysis = True #False


data_path = os.path.join(os.getcwd(), 'data')

# File paths for saving/loading data and results
fem_data_file = os.path.join(data_path, 'finite_element_model.mat')
simres_file = os.path.join(data_path, 'SimResults.json')
freqBft_file = os.path.join(data_path, 'f_stoch.mat')

# Load finite element model from a .mat file
# Expected variables inside: K, M, C, Ke, NNode, NElem, XNOD, ELEMNOA, ELEMNOB, ELEMDOF, CORRES, NDOF, iDOF_obs
fem_data = sio.loadmat(fem_data_file)
K = fem_data['K']
M = fem_data['M']
C = fem_data['C']
Ke = fem_data['Ke']
NNode = fem_data['NNode']
NElem = fem_data['NElem']
XNOD = fem_data['XNOD']
ELEMNOA = fem_data['ELEMNOA'] - 1
ELEMNOB = fem_data['ELEMNOB'] - 1
ELEMDOF = fem_data['ELEMDOF'] - 1
CORRES = fem_data['CORRES'].flatten() - 1
NDOF = int(fem_data['NDOF'].item())
iDOF_obs = fem_data['iDOF_obs']


# Check if the results file exists
if not os.path.isfile(simres_file):
    print('Results from previous analysis not found.')
    redo_analysis = True

if redo_analysis:
    
    # Define wind model for buffeting analysis
    windModel = {
        'rho': 1.22,       # Air density [kg/m³]
        'U': 34.66,        # Mean wind speed [m/s]
        'sigu': 4.56,      # Wind speed standard deviation
        'Iu': 4.56 / 34.66,  # Turbulence intensity
        'Lux': 50,         # Longitudinal turbulence scale [m]
        'C': 8             # Additional aerodynamic parameter
    }

    # Define aerodynamic section properties
    aeroSec = {
        'B': 30,          # Deck width [m]
        'CD': 0.4         # Drag coefficient at zero AoA
    }

    # Define structural model properties
    structuralModel = {
        'NDOF': NDOF,
        'M': M,
        'K': K,
        'C': C,
        'Ke': Ke,
        'nModes': 7,
        'xi_s': np.full(7, 0.003),  # Modal damping ratios - should be consistent with C
        
        'XNOD': XNOD,
        'ELEMNOA': ELEMNOA,
        'ELEMNOB': ELEMNOB,
        'ELEMDOF': ELEMDOF,
        
        'iDOF_obs': iDOF_obs,
        'corres': CORRES,
        'loadedDOFs': np.where(CORRES % 3 != 2)[0]  # Identify loaded DOFs
    }

    do_buffeting_analysis(windModel, structuralModel, aeroSec, simres_file, freqBft_file)

