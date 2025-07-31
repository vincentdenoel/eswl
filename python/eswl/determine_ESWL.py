import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh, eig
from numpy import trapz
from scipy.interpolate import interp1d
import scipy.io as sio
from eswl.tools import newmark



def determine_ESWL(time, dp, nodal_info, modal_info):
    """
    Calculate Equivalent Static Wind Loads (ESWL) for a structure based on modal analysis.
    
    Args:
        time (array): Time steps of the simulation.
        dp (array): [n x nb] Pressure time history, where n is the number of time steps, and nb is the number of pressure taps.
        nodal_info (dict): Contains structural response information for nodal loads.
            z: [nl x nr] Responses for unit loads at loading points.
        modal_info (dict): Contains modal response information.
            z: [nm x nr] Responses for unit modal amplitude.
            freq: [nm x 1] Natural frequencies of modes.
            mgen, kgen: [nm x 1] Modal masses and stiffnesses.
            damp: [nm x 1] Modal damping ratios.
            mode: [nl x 3 x nm] Mode shapes at each node.
            mil: [nl x 3 x nm] Modal inertial loads at each node.

    Returns:
        dict: A dictionary containing ESWLs, mean responses, and other structural response metrics.
    """

    # Initialize variables
    n = dp.shape[0]  # Number of time steps
    nl = nodal_info["z"].shape[0]  # Number of loading points
    nr = nodal_info["z"].shape[1]  # Number of response metrics
    nm = modal_info["z"].shape[0]  # Number of modes

    # Compute forces at nodes
    dp_nodes = dp.T  # Transpose to get forces at nodes [nb x n]

    # Separate mean and fluctuation components
    mean_dp = np.mean(dp_nodes, axis=1)  # Mean pressure
    mean_F = mean_dp  # Mean force at nodes
    fluct_dp = dp_nodes - mean_dp[:, None]  # Fluctuations in pressure
    fluct_F = fluct_dp  # Fluctuations in force

    # Covariance of forces at nodes
    cov_F = np.cov(fluct_F)  # Covariance of fluctuating forces
    var_F = np.diag(cov_F)  # Variances of forces
    cov_F_coh = np.sqrt(np.outer(var_F, var_F))  # Fully coherent covariance matrix

    # Background responses (nodal basis)
    zn = nodal_info["z"]  # Static responses for unit loads
    mean_z = zn.T @ mean_F  # Mean responses
    zb = zn.T @ fluct_F  # Background components of responses
    cov_z = zn.T @ cov_F @ zn  # Covariance of background responses
    cov_z_coh = zn.T @ cov_F_coh @ zn  # Fully coherent background covariance

    # Resonant responses (modal basis)
    zm = modal_info["z"]  # Modal responses for unit amplitudes
    F_star = np.zeros((nm, n))  # Modal excitation forces

    for mode in range(nm):
        for node in range(nl):
            F_star[mode, :] += fluct_F[node, :] * modal_info["mode"][node, 1, mode]

    dt = time[1] - time[0]  # Time step size
    qb = np.zeros((nm, n))
    q = np.zeros((nm, n))
    qr = np.zeros((nm, n))
    zr_mode = np.zeros((n, nr, nm))
    std_zr_mode = np.zeros((nr, nm))

    for mode in range(nm):
        f = modal_info["freq"][mode]
        m = modal_info["mgen"][mode]
        xi = modal_info["damp"][mode]
        p = F_star[mode, :]
        time, dep, _, _, dep_st = newmark(m, xi, f, p, dt, n)
        qb[mode, :] = dep_st[0:n]  # Background modal response
        q[mode, :] = dep[0:n]  # Total modal response
        qr[mode, :] = dep[0:n] - dep_st[0:n]  # Resonant part of modal response
        zr_mode[:, :, mode] = qr[mode, 0:n, None] @ zm[mode, None, 0:n]
        std_zr_mode[:, mode] = np.std(zr_mode[:, :, mode], axis=0)

    zr = np.sum(zr_mode, axis=2).T  # Total resonant response
    z = zr + zb  # Total response (background + resonant)

    # Standard deviations
    std_q = np.std(qr, axis=1)
    std_z = np.std(z, axis=1)
    std_zb = np.std(zb, axis=1)
    std_zr = np.std(zr, axis=1)

    # Extreme value statistics
    T = 600  # Duration for extreme value analysis
    nw = int(T / dt)  # Number of windows
    nB = max(1, n // nw)

    zminBlock = np.zeros((nr, nB))
    zmaxBlock = np.zeros((nr, nB))

    for block in range(nB):
        start = block * nw
        end = (block + 1) * nw
        zminBlock[:, block] = np.min(z[:, start:end], axis=1)
        zmaxBlock[:, block] = np.max(z[:, start:end], axis=1)

    zMin_mean = np.mean(zminBlock, axis=1)
    zMax_mean = np.mean(zmaxBlock, axis=1)

    gmin = zMin_mean / std_z
    gmax = zMax_mean / std_z
    g = (gmax - gmin) / 2
    z_max = mean_z + g * std_z
    z_min = mean_z - g * std_z

    # Compute ESWL
    eswl_b = np.zeros((nl, nr))
    eswl_b_coh = np.zeros((nl, nr))
    ESWL_b = np.zeros((nl, nr, 3))
    ESWL_r = np.zeros((nl, nr, 3))
    ESWL = np.zeros((nl, nr, 3))
    wght_b = np.zeros(nr)
    wght_r = np.zeros((nm, nr))
    eswl_r = np.zeros((nl, nr, 3))
    z_rec = np.zeros((nr, nr))

    for r in range(nr):
        # Magnitude of force (normal to surface)
        eswl_b[:, r] = g[r] * nodal_info['z'][:, r].T @ cov_F / std_zb[r]
        eswl_b_coh[:, r] = g[r] * nodal_info['z'][:, r].T @ cov_F_coh / std_zb[r]

        wght_b[r] = std_zb[r] / std_z[r]

        # Components of the force in each direction (xyz)
        ESWL_b[:, r] = eswl_b[:, r][:, None] * wght_b[r]  # weighting

        # Reconstructed responses (without average responses, in mean_z, to be added)
        z_rec[:, r] = nodal_info['z'].T @ (eswl_b[:, r] * wght_b[r])  # background

        for mode in range(nm):
            wght_r[mode, r] = std_zr_mode[r, mode] / std_z[r]

            # Ensure dynamic contributions are positive
            wght_r[mode, r] *= np.sign(modal_info['z'][mode, r])

            eswl_r[:, r, 0] = g[r] * std_q[mode] * modal_info['mil'][:, 0, mode]  # along x
            eswl_r[:, r, 1] = g[r] * std_q[mode] * modal_info['mil'][:, 1, mode]  # along y
            eswl_r[:, r, 2] = g[r] * std_q[mode] * modal_info['mil'][:, 2, mode]  # along z

            ESWL_r[:, r] += eswl_r[:, r] * wght_r[mode, r]  # weighting

            z_rec[:, r] += modal_info['z'][mode, :].T * (g[r] * std_q[mode] * wght_r[mode, r])  # resonant

    # Final ESWL calculation
    ESWL = ESWL_b + ESWL_r

    # Mean forces
    ESWL_m = mean_F

    return {
        "mean_z": mean_z,
        "std_z": std_z,
        "std_zb": std_zb,
        "std_zr": std_zr,
        "z": z,
        "z_min": z_min,
        "z_max": z_max,
        "z_rec": z_rec,
        "ESWL_m": ESWL_m,
        "ESWL_b": ESWL_b,
        "ESWL_r": ESWL_r,
        "ESWL": ESWL,
    }

