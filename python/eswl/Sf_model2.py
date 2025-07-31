import numpy as np


def Sf_model2(omega, x, windModel, aeroSec):
    """
    Returns the x-PSD matrix of drag loads.
    
    Parameters:
    - omega: array of angular frequencies (rad/s), shape (n,)
    - x: array of positions along the structure, shape (nx,)
    - windModel: object or dict with fields:
        * U   : mean wind velocity (m/s)
        * Iu  : turbulence intensity (dimensionless)
        * Lux : turbulence lengthscale (m)
        * rho : air density (kg/m^3)
        * C   : coherence decay constant (dimensionless)
    - aeroSec: object or dict with fields:
        * B  : characteristic width of section (m)
        * CD : drag coefficient (dimensionless)
        
    Returns:
    - Sf_drag: PSD matrix of drag loads, shape (nx, nx, n)
    """
    
    nx = len(x)
    omega = np.atleast_1d(omega)
    n = len(omega)

    Sf_drag = np.zeros((nx, nx, n))

    rho = windModel['rho']    # kg/m3
    U = windModel['U']        # m/s
    Iu = windModel['Iu']
    L = windModel['Lux']
    sig_u = Iu * U

    B = aeroSec['B']          # m

    q = 0.5 * rho * U**2      # N/m^2
    C = windModel['C']        # dimensionless

    # Matrix of distances between nodes (nx x nx)
    dx = np.abs(np.subtract.outer(x, x))

    # Frequency (Hz)
    f = omega / (2 * np.pi)

    # Monin frequency
    y = f * L / U

    # PSD of turbulence velocity component u(t)
    Su = 4 * L / U * sig_u**2 / (1 + 70.7 * y**2)**(5/6)  # shape (n,)

    # Aerodynamic coefficient
    coef_D_u = q * B * (2 * aeroSec['CD'] / U)  # Ns/m^2

    # Consider exponential coherence
    for i in range(n):
        coh = np.exp(-C * omega[i] * dx / (2 * np.pi * U))  # shape (nx, nx)
        Sf_drag[:, :, i] = coef_D_u**2 * Su[i] * coh          # N^2 s / m^2

    # Multiply by length of each element
    dx2 = (x[1:] - x[:-1]) / 2
    DX = np.concatenate((dx2, [0])) + np.concatenate(([0], dx2))
    DX = np.outer(DX, DX)
    
    for i in range(n):
        Sf_drag[:, :, i] *= DX  # N^2 s

    return Sf_drag
