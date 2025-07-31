import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh, eig
from numpy import trapz
from scipy.interpolate import interp1d
import scipy.io as sio


def gen_samples(spectral_density_handle, time, eigenfreqs):
    """
    Generates samples from a given spectral density matrix.

    Parameters:
    spectral_density_handle: Callable
        Function handle to compute spectral density matrix.
    time: numpy.ndarray
        Time vector.
    eigenfreqs: numpy.ndarray
        Frequencies for eigenvalue decomposition.

    Returns:
    time: numpy.ndarray
        Modified time vector.
    samples: numpy.ndarray
        Generated samples in the time domain.
    """
    tmp = spectral_density_handle(0).squeeze()
    ndof = np.size(tmp,1)
    cov_smp = np.zeros((ndof, ndof), dtype=float)

    # 1. Time vector adjustment for even length
    if len(time) % 2:
        time = np.append(time, 2 * time[-1] - time[-2])
    n = len(time)
    dt = time[1] - time[0]
    df = 1 / (time[-1])

    time = np.arange(0, n) * dt
    freqs = np.arange(0, n) * df

    if eigenfreqs[-1] < freqs[-1]:
        eigenfreqs = np.append(eigenfreqs, freqs[-1])

    # 2. Modal decomposition of the spectral density matrix
    modal_matrix = np.zeros((ndof, ndof, len(eigenfreqs)), dtype=float)
    eigenvalue_matrix = np.zeros((ndof, len(eigenfreqs)), dtype=float)

    for freq_idx in range(len(eigenfreqs)):
        spec_density = spectral_density_handle(eigenfreqs[freq_idx] * 2 * np.pi).squeeze()        
        eigvals, eigvecs = np.linalg.eigh(spec_density)

        sorted_indices = np.argsort(eigvals)[::-1]
        eigvals = eigvals[sorted_indices]
        eigvecs = eigvecs[:, sorted_indices]

        if freq_idx > 0:
            cov_smp += spec_density * (eigenfreqs[freq_idx] - eigenfreqs[freq_idx - 1]) * (4 * np.pi)
            for mode_idx in range(eigvecs.shape[1]):
                if np.dot(prev_eigvecs[:, mode_idx], eigvecs[:, mode_idx]) < 0:
                    eigvecs[:, mode_idx] *= -1                
        else:
            cov_smp += spec_density * (eigenfreqs[1] - eigenfreqs[0]) * (4 * np.pi)

        prev_eigvecs = eigvecs
        modal_matrix[:, :, freq_idx] = eigvecs
        eigenvalue_matrix[:, freq_idx] = eigvals    
    
    #cov_modal = trapz(eigenvalue_matrix * 4 * np.pi, eigenfreqs)

    # 3. Generate modal processes
    modal_process = np.zeros((n, ndof), dtype=complex)

    for mode_idx in range(ndof):
        interpolated_eigenvalues = interp1d(eigenfreqs, eigenvalue_matrix[mode_idx, :], kind='linear', fill_value="extrapolate")(freqs)
        random_phases = np.random.uniform(0, 2 * np.pi, n // 2 - 1)
        scaling_factor = 4 * np.pi / (2 * dt)

        modal_process[0, mode_idx] = 0
        modal_process[1:n // 2, mode_idx] = scaling_factor * np.sqrt(n * dt / (2 * np.pi) * interpolated_eigenvalues[1:n // 2]) * \
                                             np.exp(1j * random_phases)
        modal_process[n // 2, mode_idx] = scaling_factor * np.sqrt(n * dt / (2 * np.pi) * interpolated_eigenvalues[n // 2])
        modal_process[n // 2 + 1:, mode_idx] = np.conj(modal_process[n // 2 - 1:0:-1, mode_idx])

    # modal_samples = np.fft.ifft(modal_process, axis=0)

    # 4. Generate nodal processes
    nodal_process = np.zeros((n, ndof), dtype=complex)


    for idof in range(ndof):
        nodal_matrix = np.zeros((n, ndof), dtype=complex)

        for mode_idx in range(ndof):
            i = np.arange(1, n//2)

            interp_func = interp1d(eigenfreqs, modal_matrix[idof, mode_idx, :], kind='linear', bounds_error=False, fill_value=0)
            nodal_matrix[i, mode_idx] = interp_func(freqs[i])
            nodal_matrix[n//2, mode_idx] = np.real(interp_func(freqs[n//2]))
            nodal_matrix[n//2 + 1:, mode_idx] = np.conj(nodal_matrix[n//2 - 1:0:-1, mode_idx])

        nodal_process[:, idof] = np.sum(nodal_matrix * modal_process, axis=1)

    # 5. Perform inverse FFT to get time-domain samples
    samples = np.fft.ifft(nodal_process, axis=0).real

    return time, samples


# Helper to test vector input support
def is_vector_input_supported(func_handle):
    """
    Check if the function handle supports vector input.

    Parameters:
    func_handle: Callable
        Function handle to test.

    Returns:
    bool: True if vector input is supported, False otherwise.
    """
    try:
        func_handle(np.array([1, 2, 3, 4]))
        return True
    except Exception:
        return False

