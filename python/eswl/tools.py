import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.io as sio


def internal_forces(NDOF, X, Ke, CORRES, NElem, ELEMDOF):
    """
    Computes the internal forces (Axial, Shear, Bending) for each element.

    Parameters:
    - NDOF: Total number of degrees of freedom (3 dofs per node, only 2 are used in this example)
    - X: Displacement vector
    - Ke: Element stiffness matrices (3D array)
    - CORRES: Mapping of global to local degrees of freedom
    - NElem: Number of elements
    - ELEMDOF: Degrees of freedom for each element (2D array)

    Returns:
    - Axial: Axial forces for each element
    - Shear: Shear forces for each element
    - Bendi: Bending moments for each element
    """
    XX = np.zeros(NDOF)
    XX[CORRES] = X.ravel()

    Axial = np.zeros((NElem, 2))
    Shear = np.zeros((NElem, 2))
    Bendi = np.zeros((NElem, 2))

    for iel in range(NElem):
        xloc = XX[ELEMDOF[:, iel]]
        Fint = Ke[:, :, iel] @ xloc
        Axial[iel, 0] = Fint[0]
        Axial[iel, 1] = -Fint[3]
        Shear[iel, 0] = Fint[1]
        Shear[iel, 1] = -Fint[4]
        Bendi[iel, 0] = -Fint[2]
        Bendi[iel, 1] = Fint[5]

    return Axial, Shear, Bendi


def FRF_mdof(M, K, C, f):
    f = np.atleast_1d(f)
    Nfr = len(f)
    N = np.size(M,0)
    H = np.zeros((N, N, Nfr), dtype=complex)    
    for i in range(Nfr):
        omega = 2 * np.pi * f[i]
        H[:, :, i] = np.linalg.inv(-M * omega**2 + 1j * omega * C + K)
    return H

def NormalizeMode(V, normMode=1):
    """
    Normalizes the matrix of modes V to a maximum real unitary for each mode.

    Parameters:
    - V: 2D numpy array of modes
    - normMode: Normalization mode (1 for unitary normalization, else for norm-based)

    Returns:
    - V_normalized: Normalized mode matrix
    """
    N = V.shape[1]
    V_normalized = np.zeros_like(V, dtype=complex)

    for s in range(N):
        vec = V[:, s]
        veca = np.abs(vec)**2  # Element-wise square of absolute values
        imax = np.argmax(veca)  # Index of maximum value

        if normMode == 1:
            vec = vec / vec[imax]
        else:
            vec = vec / np.linalg.norm(vec)
        
        V_normalized[:, s] = vec

    return V_normalized



def plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, eswl=None, fit_env_scale=False, doplot=True, ax=None):
    """
    Plots displacement, bending moment, and static wind loading diagrams.

    Parameters:
    - ESWLresults: Results object with fields `z_max`, `z_min`, `ESWL`, `z_rec`
    - structuralModel: Structural model containing XNOD, K, Ke, corres, ELEMNOA, ELEMDOF
    - i_x: Indices in z with displacements and rotations of nodes
    - i_m: Indices in z with bending moments
    - i_z: List of responses
    - clr: Color for plotting
    - eswl: Optional ESWL matrix
    - fit_env_scale: Whether to fit envelope scale (default False)
    - doplot: Whether to produce plots (default True)
    """
    XNOD = structuralModel['XNOD']
    upper = ESWLresults['z_max']
    lower = ESWLresults['z_min']

    if eswl is not None:
        if eswl.ndim == 1:
            eswl = eswl[:, np.newaxis]
            
        # Recompute displacement and bending moment for a given ESWL and scale it
        K = structuralModel['K']

        f = np.zeros((eswl.shape[0] * 2, len(i_z)))
        f[0::2, :] = eswl[:, i_z]
        xx = np.linalg.solve(K, f)

        _, _, Bendi = internal_forces(
            3*len(XNOD), xx, structuralModel['Ke'],
            structuralModel['corres'], len(structuralModel['ELEMNOA'][0]),
            structuralModel['ELEMDOF']
        )

        if fit_env_scale:
            scale_x = np.min(np.abs(upper[i_x] / xx[0::2, 0]))
            scale_m = np.min(np.abs(upper[i_m[:-1]] / Bendi[:, 0]))
            scale = min(scale_x, scale_m)
        else:
            scale = upper[i_z] / xx[i_z]

        xx *= scale
        Bendi *= scale
        eswl *= scale

        x_rec = xx[0::2, :]
        m_rec = np.concatenate([Bendi[:, 0], [Bendi[-1, 1]]])

    else:
        eswl = ESWLresults['ESWL'][:, :, 2]
        x_rec = ESWLresults['z_rec'][i_x, i_z]
        m_rec = ESWLresults['z_rec'][i_m, i_z]
        scale = 1

    if not doplot:
        return x_rec, m_rec, scale

    if np.any(np.array(i_z) <= 170):
        i__z = (np.array(i_z) + 1) // 2  # Displacement
    else:
        i__z = np.array(i_z) - 170  # Bending moment

    plotSym = True

    if ax is None:
        fig, ax = plt.subplots(3, 1, figsize=(10, 6))

    # Plot Displacements
    ax[0].plot(XNOD, upper[i_x], 'k:', linewidth=1.2)
    ax[0].plot(XNOD, lower[i_x], 'k:', linewidth=1.2)
    ax[0].plot(XNOD, x_rec, color=clr)
    if plotSym:
        ax[0].plot(XNOD, -x_rec, color=clr)
    ax[0].plot(XNOD[i__z], x_rec[i__z], 'ko', markerfacecolor='k')
    ax[0].set_title('Displacements')
    ax[0].set_xlim(XNOD[0], XNOD[-1])
    ax[0].axis('off')

    # Plot Bending Moments
    ax[1].plot(XNOD, upper[i_m], 'k:', linewidth=1.2)
    ax[1].plot(XNOD, lower[i_m], 'k:', linewidth=1.2)
    ax[1].plot(XNOD, m_rec, color=clr)
    if plotSym:
        ax[1].plot(XNOD, -m_rec, color=clr)
    ax[1].plot(XNOD[i__z], m_rec[i__z], 'ko', markerfacecolor='k')
    ax[1].set_title('Bending Moment [kNm]')
    ax[1].set_xlim(XNOD[0], XNOD[-1])
    ax[1].axis('off')

    # Plot Static Wind Loading
    ax[2].plot(XNOD, eswl[:, i_z], color=clr)
    ax[2].plot([XNOD[0], XNOD[-1]], [0, 0], 'k', linewidth=0.3)
    ax[2].set_title('Static wind loading')
    ax[2].set_xlim(XNOD[0], XNOD[-1])
    ax[2].axis('off')

    plt.tight_layout()

    return x_rec, m_rec, scale


def plotTF(t, x, NFFT=512, sp=None, scaling=1, *varargin):
    """
    Plots time series and PSD for given signal data.
    
    Parameters:
    - t: Time array
    - x: Signal array
    - NFFT: Number of FFT points (default 512)
    - sp: Subplot axes (default creates two subplots)
    - scaling: Scaling factor for PSD (default 1)
    - varargin: Additional plot options for the time series
    
    Returns:
    - f: Frequency array
    - psdX: PSD matrix of the signals
    """
    if x.shape[0] > x.shape[1]:
        x = x.T  # Ensure x is in shape (signal, samples)

    dt = t[1] - t[0]

    XX, f, _ = VincePSD(x, NFFT, dt, scaling)

    psdX = np.zeros((XX.shape[0], XX.shape[1]))
    for i in range(x.shape[0]):
        psdX[:, i] = XX[:, i, i]
        psdX[0, i] = np.nan

    # Create subplots if not provided
    if sp is None:
        fig, sp = plt.subplots(2, 1)

    # Plot time series
    ax1 = sp[0]
    ax1.plot(t, x.T, *varargin)
    ax1.set_xlabel('Time [s]')
    ax1.set_title('Time Series')

    # Plot PSD
    ax2 = sp[1]
    ax2.semilogy(f, psdX)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_title('Power Spectral Density')

    plt.tight_layout()
    plt.show()

    return f, psdX

def VincePSD(x, NFFT, DT, scaling=1):
    """
    Computes the Power Spectral Density (PSD) matrix of signals.
    
    Parameters:
    - x: Input signal array (shape: signals x samples)
    - NFFT: Number of FFT points
    - DT: Time step
    - scaling: Scaling factor for PSD (default 1)
    
    Returns:
    - PSD: PSD matrix
    - FREQ: Frequency array
    - OMEGA: Angular frequency array
    """
    NS, N = x.shape

    # Replace NaNs in signals with the mean of the signal
    for i in range(NS):
        x[i, np.isnan(x[i, :])] = np.nanmean(x[i, :])

    Nblocs = N // NFFT
    PSD = np.zeros((NFFT // 2, NS, NS))

    # Create Hanning window
    hann = np.sin(np.pi * np.arange(NFFT) / (NFFT - 1)) ** 2
    W = np.sum(hann ** 2)

    # Loop over blocks
    for i in range(Nblocs):
        for s1 in range(NS):
            xx1 = x[s1, i * NFFT:(i + 1) * NFFT] * hann
            XX1 = np.fft.fft(xx1 - np.mean(xx1))
            for s2 in range(NS):
                xx2 = x[s2, i * NFFT:(i + 1) * NFFT] * hann
                XX2 = np.fft.fft(xx2 - np.mean(xx2))
                periodogram = XX1 * np.conj(XX2)
                PSD[:, s1, s2] += periodogram[:NFFT // 2].real

    PSD /= Nblocs * W

    DF = 1 / (NFFT * DT)
    FREQ = np.arange(NFFT // 2) * DF
    OMEGA = 2 * np.pi * FREQ

    if scaling == 1:
        PSD = PSD / NFFT / DF * 2
    elif scaling == 2:
        PSD = PSD / NFFT / DF * 2 / (4 * np.pi)

    return PSD, FREQ, OMEGA

def newmark(m, xi, f, p, dt, Nstep, a=0.25, d=0.5, dep0=0, vit0=0):
    
    omega = 2 * np.pi * f
    k = m * omega ** 2
    c = 2 * m * omega * xi

    dep = np.zeros(Nstep+1)
    dep[0] = dep0
    vit = np.zeros(Nstep+1)
    vit[0] = vit0
    acc = np.zeros(Nstep+1)
    qs = np.zeros(Nstep+1)
    akf = m / a / dt ** 2 + d / a / dt * c + k

    for i in range(1, Nstep):
        dep[i] = (p[i] + m * (1 / a / dt ** 2 * dep[i-1] + 1 / a / dt * vit[i-1] + (1 / 2 / a - 1) * acc[i-1])
                + c * (d / a / dt * dep[i-1] + (d / a - 1) * vit[i-1] + dt / 2 * (d / a - 2) * acc[i-1])) / akf
        vit[i] = d / a / dt * (dep[i] - dep[i-1]) + (1 - d / a) * vit[i-1] + dt * (1 - d / 2 / a) * acc[i-1]
        acc[i] = 1 / a / dt ** 2 * (dep[i] - dep[i-1]) - 1 / a / dt * vit[i-1] - (1 / 2 / a - 1) * acc[i-1]
        qs[i] = p[i] / k

    time = np.linspace(0.,Nstep*dt,Nstep+1)
    return time, dep, vit, acc, qs


def serialize(obj):
    """
    Recursively serialize objects to make them JSON-compatible.
    Converts numpy arrays to lists and handles other unsupported types.
    """
def serialize(obj):
    if isinstance(obj, np.ndarray):
        if np.iscomplexobj(obj):
            # Serialize complex ndarray as a dict with separate real and imag parts
            return {"type": "complex_ndarray", "real": obj.real.tolist(), "imag": obj.imag.tolist()}        
        return obj.tolist()  # Convert NumPy arrays to lists
    elif isinstance(obj, dict):
        return {key: serialize(value) for key, value in obj.items()}  # Recursively process dicts
    elif isinstance(obj, (list, tuple)):
        return [serialize(item) for item in obj]  # Recursively process lists/tuples
    elif isinstance(obj, (np.integer, np.floating)):  # Handle NumPy scalar types
        return obj.item()
    else:
        return obj  # Return the object itself if it's already JSON-serializable

def deserialize(obj):
    """
    Recursively deserialize objects from JSON-compatible format.
    Converts lists back to NumPy arrays where appropriate, including complex ndarrays.
    """
    if isinstance(obj, dict):
        # Check if the dict represents a complex ndarray
        if obj.get("type") == "complex_ndarray":
            real = np.array(obj["real"])
            imag = np.array(obj["imag"])
            return real + 1j * imag
        # Check if the dict represents a single complex number
        if obj.get("type") == "complex":
            return complex(obj["real"], obj["imag"])
        # Recursively process other dicts
        return {key: deserialize(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return np.array(obj)  # Convert regular lists back to NumPy arrays
    else:
        return obj  # Return the object itself if it's not a list or dict
