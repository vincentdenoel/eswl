import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh, eig
from numpy import trapz
from scipy.interpolate import interp1d
import scipy.io as sio
import json

from eswl.tools import internal_forces, FRF_mdof
from eswl.tools import serialize
from eswl.Sf_model2 import Sf_model2
from eswl.gen_samples import gen_samples
from eswl.determine_ESWL import determine_ESWL

def do_buffeting_analysis(windModel, structuralModel, aeroSec, simres_file, freqBft_file):
    """
    Perform buffeting analysis of a finite element model. 
    This script sequentially:
    - Runs spectral analysis in the frequency domain.
    - Runs spectral analysis in the time domain (for illustration purposes).
    - Determines ESWLs for displacements and bending moments.

    Args:
    - windModel: Turbulence model (used for buffeting analysis).
    - structuralModel: Contains mode shapes, damping, modal masses, etc.
    - aeroSec: Static aerodynamic coefficients (Drag only is used in this example).
    - simres_file: String containing the file name to save results.
    """

    results = {}

    # Structural properties (mass, stiffness and damping matrices / damping ratio)
    M = structuralModel['M']
    K = structuralModel['K']
    C = structuralModel['C']
    xi_s = structuralModel['xi_s']

    # Node coordinates (1-D in this example), number of DOFs, number of modes for modal analysis 
    XNOD = structuralModel['XNOD']
    NDOF = structuralModel['NDOF']
    nModes = structuralModel['nModes']
    NElem = len(structuralModel['ELEMNOA'][0])

    # Compute mode shapes
    eigvals, eigvecs = eigh(K, M)
    V = eigvecs[:, :nModes]
    V = V / np.max(np.abs(V), axis=0)    
    omega = np.sqrt(eigvals[:nModes])
    fnat = omega / (2 * np.pi)

    Mgen = np.diag(V.T @ M @ V)
    Kgen = np.diag(V.T @ K @ V)
    Cgen = 2 * np.sqrt(Mgen * Kgen) * xi_s

    A_disp = V  # Influence matrix for displacements (mode shapes)

    ABendi = []
    for imod in range(nModes):
        _, _, Bendi = internal_forces(3*len(XNOD), V[:, imod], structuralModel['Ke'],
                                    structuralModel['corres'], NElem, structuralModel['ELEMDOF'])
        ABendi.append(np.hstack((Bendi[:, 0], Bendi[-1, 1])))
    ABendi = np.array(ABendi).T

    #ABendi = ABendi[0:-1, :] # remove last entry (it is zero)

    A_Bendi = []
    for idof in range(len(K)):
        x_tmp = np.zeros(len(K))
        x_tmp[idof] = 1
        _, _, Bendi = internal_forces(3*len(XNOD), x_tmp, structuralModel['Ke'],
                                    structuralModel['corres'], NElem, structuralModel['ELEMDOF'])
        A_Bendi.append(np.hstack((Bendi[:, 0], Bendi[-1, 1])))
    A_Bendi = np.array(A_Bendi).T

    # 1. Spectral analysis (frequency domain)
    print(' * Running spectral analysis in frequency domain...')    
    f_stoch = sio.loadmat(freqBft_file)['f_stoch'][0]

    v = V[structuralModel['loadedDOFs'], :]
    kn = np.arange(v.shape[0])

    # initialize empty matrices to store PSDs of modal and nodal responses
    Sq = np.zeros((nModes, nModes, len(f_stoch)), dtype=complex)
    SqB = np.zeros((nModes, nModes, len(f_stoch)), dtype=complex)
    Sx = np.zeros((NDOF, NDOF, len(f_stoch)), dtype=complex)
    Sxf = np.zeros((NDOF, len(kn), len(f_stoch)), dtype=complex)
    Sf = np.zeros((len(kn), len(kn), len(f_stoch)), dtype=complex)

    for ifr, f in enumerate(f_stoch):
        w = 2 * np.pi * f
        x = XNOD.flatten()
        Sf_drag = Sf_model2(w, x, windModel, aeroSec).squeeze()
        Sfgen = v.T @ Sf_drag[np.ix_(kn, kn)] @ v

        H = FRF_mdof(np.diag(Mgen), np.diag(Kgen), np.diag(Cgen), f).squeeze()
        Sq[:, :, ifr] = H @ Sfgen @ H.conj().T
        SqB[:, :, ifr] = np.outer(1 / Kgen, np.diag(Sfgen) / Kgen)

        H__ = FRF_mdof(M, K, C, f)
        H = H__[:, structuralModel['loadedDOFs']].squeeze()
        Sx[:, :, ifr] = H @ Sf_drag[np.ix_(kn, kn)] @ H.conj().T
        Sxf[:, :, ifr] = H @ Sf_drag[np.ix_(kn, kn)]
        Sf[:, :, ifr] = Sf_drag

    # Modal basis responses
    covq = trapz(Sq, f_stoch, axis=2).real
    covqB = trapz(SqB, f_stoch, axis=2).real
    covx = A_disp @ covq @ A_disp.T
    covf = trapz(Sf, f_stoch, axis=2).real
    stdx = np.sqrt(np.diag(covx))
    covxf = trapz(Sxf, f_stoch, axis=2).real
    covbending = ABendi @ covq @ ABendi.T
    stdbending = np.sqrt(np.diag(covbending))

    results['nodal'] = {
        'covq': covq,
        'cov_x': covx,
        'cov_f': covf,
        'cov_xf': covxf,
        'std_x': stdx,
        'cov_M': covbending,
        'std_M': stdbending
    }

    # 2. Buffeting analysis in time domain
    print(' * Running buffeting analysis in time domain...')
    dragforces_file = os.path.join(os.getcwd(), 'data', 'drag_forces.json')
    if not os.path.exists(dragforces_file):
        print('Drag forces not found. Generating...')
        dt = 0.04
        S_hdl = lambda w: Sf_model2(w, x, windModel, aeroSec) / (4 * np.pi)
        time, drag_forces = gen_samples(S_hdl, np.arange(2**16) * dt, np.concatenate(([0], 2**np.linspace(-9, 4, 261))))
        data = {
            'time': time.flatten().tolist(),
            'drag_forces': drag_forces.tolist()
        }
        print('Done. Saving to json file for future use...')
        with open(dragforces_file, 'w') as json_file:
            json.dump(data, json_file, indent=4)
        print('Done.')
    else:
        print('Reloading time series of wind loads...')
        with open(dragforces_file, 'r') as json_file:
            data = json.load(json_file)
        time = np.array(data['time'])
        drag_forces = np.array(data['drag_forces'])
        print('Done.')
        
    nl = drag_forces.shape[1]  # Number of loading points

    # Temporary plot of drag forces, uncomment if needed
    #from eswl.tools import plotTF
    #plotTF(time, drag_forces[:,0:2], NFFT=2048*4)


    # Initialize nodalInfo
    nodalInfo = {
        'f_stoch': f_stoch,
        'Sx': Sx,
        'z_x': [],
        'z_m': [],
    }

    # Responses under unit loads at each 'pressure tap'
    for i in range(nl):
        p = np.zeros(K.shape[0])
        p[2 * i] = 1
        x = np.linalg.solve(K, p)  # Solve Kx = p

        Axial, Shear, Bendi = internal_forces(
            3*len(XNOD), x, structuralModel['Ke'], structuralModel['corres'],
            NElem, structuralModel['ELEMDOF']
        )
        
        nodalInfo['z_x'].append(x)
        nodalInfo['z_m'].append(np.concatenate([Bendi[:, 0]]))

    nodalInfo['z_x'] = np.array(nodalInfo['z_x'])
    nodalInfo['z_m'] =  np.hstack([np.array(nodalInfo['z_m']), np.zeros((np.array(nodalInfo['z_m']).shape[0], 1))])
    nodalInfo['z'] = np.hstack([nodalInfo['z_x'], nodalInfo['z_m']])

    # Modal information
    mil = K @ V
    modalInfo = {
        'z_x': V.T,
        'z_m': ABendi.T,
        'z': np.hstack([V.T, ABendi.T]),
        'freq': fnat,
        'mgen': Mgen,
        'kgen': Kgen,
        'nm': nModes,
        'damp': xi_s,
        'mode': np.zeros((V.shape[0] // 2, 3, V.shape[1])),
        'mil': np.zeros((mil.shape[0] // 2, 3, mil.shape[1])),
    }

    # Populate modal info, modes
    modalInfo['mode'][:, 1, :] = V[::2, :]
    modalInfo['mil'][:, 1, :] = mil[::2, :]

    # 3. Compute ESWL results
    print ('Determining ESWLs...')
    ESWLresults = determine_ESWL(time, drag_forces, nodalInfo, modalInfo)

    # Save results
    print ('Saving results to', simres_file, '...')
    data_to_save = {
    'results': results,
    'ESWLresults': ESWLresults,
    'modalInfo': modalInfo,
    'nodalInfo': nodalInfo,
    'time': time,
    'drag_forces': drag_forces,
    'structuralModel': structuralModel,
    'windModel': windModel,
    'aeroSec': aeroSec,
    }
    with open(simres_file, 'w') as json_file:
        json.dump(serialize(data_to_save), json_file, indent=4)

    print('Done.')
