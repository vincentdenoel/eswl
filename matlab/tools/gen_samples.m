function [time, samples] = gen_samples(spectral_density_handle, time, eigenfreqs)
% GEN_SAMPLES generates samples from a given spectral density matrix
%
% [time, samples] = gen_samples(spectral_density_handle, time, eigenfreqs)
% Generates time-domain samples based on the input spectral density.
% 
% spectral_density_handle: Function handle to spectral density matrix
% time: Time vector
% eigenfreqs: Frequencies for which eigen value decomposition is performed
%
% Output:
% time: Modified time vector
% samples: Generated samples
%

%% 1. Time vector adjustment for even length
if mod(length(time), 2)
    time(end+1) = 2 * time(end) - time(end-1);
end
n = length(time);

df = 1 / time(end);
dt = time(2) - time(1);

time = (0:n-1) * dt; % Recalculate time vector
freqs = (0:n-1) * df; % Frequency vector

if eigenfreqs(end) < freqs(end)
    eigenfreqs(end+1) = freqs(end); % Extend frequency range
end

%% 2. Modal decomposition of the spectral density matrix
% Compute modal decomposition for each frequency and interpolate if needed
if isVectorInputSupported(spectral_density_handle)
    spectral_density = spectral_density_handle(eigenfreqs * 2 * pi); % Use vector input

    progress_bar = waitbar(0, 'Performing Cholesky decomposition...');
    for freq_idx = 1:length(eigenfreqs)
        [eigenvectors, eigenvalues] = eig(spectral_density(:, :, freq_idx));
        eigenvalues_diag = diag(eigenvalues);
        [eigenvalues_diag, order] = sort(eigenvalues_diag, 'descend');
        eigenvectors = eigenvectors(:, order);

        % Ensure consistency of eigenvector orientation
        if freq_idx > 1
            for mode_idx = 1:size(eigenvectors, 2)
                if eigenvectors(:, mode_idx)' * prev_eigenvectors(:, mode_idx) < 0
                    eigenvectors(:, mode_idx) = -eigenvectors(:, mode_idx);
                end
            end
        end
        prev_eigenvectors = eigenvectors;
        modal_matrix(:, :, freq_idx) = eigenvectors;
        eigenvalue_matrix(:, freq_idx) = eigenvalues_diag;
        waitbar(freq_idx / length(eigenfreqs), progress_bar);
    end
    delete(progress_bar);
    
else
    progress_bar = waitbar(0, 'Performing Cholesky decomposition...');
    for freq_idx = 1:length(eigenfreqs)
        spectral_density = spectral_density_handle(eigenfreqs(freq_idx) * 2 * pi);
        [eigenvectors, eigenvalues] = eig(spectral_density);
        eigenvalues_diag = diag(eigenvalues);
        [eigenvalues_diag, order] = sort(eigenvalues_diag, 'descend');
        eigenvectors = eigenvectors(:, order);

        modal_matrix(:, :, freq_idx) = eigenvectors;
        eigenvalue_matrix(:, freq_idx) = eigenvalues_diag;
        waitbar(freq_idx / length(eigenfreqs), progress_bar);
    end
    delete(progress_bar);
end

target_covariance = trapz(4 * pi * eigenfreqs, spectral_density, 3);

num_histories = size(spectral_density, 1); % Number of processes
num_modes = num_histories; % Number of modes

%% 3. Generate modal processes
modal_process = zeros(n, num_modes);

for mode_idx = 1:num_modes
    interpolated_eigenvalues = interp1(eigenfreqs, eigenvalue_matrix(mode_idx, :), freqs);
    random_phases = 2 * pi * rand(1, n / 2 - 1);
    scaling_factor = 4 * pi / (2 * dt);

    modal_process(1, mode_idx) = 0;
    i = 2:n/2;
    modal_process(i, mode_idx) = scaling_factor * sqrt(n * dt / (2 * pi) * interpolated_eigenvalues(i)) .* exp(1i * random_phases);

    modal_process(n/2+1, mode_idx) = scaling_factor * sqrt(n * dt / (2 * pi) * interpolated_eigenvalues(n/2+1));
    modal_process(n/2+2:end, mode_idx) = conj(modal_process(n/2:-1:2, mode_idx));
end

%% 4. Generate nodal processes
nodal_matrix = zeros(n, num_modes);

for history_idx = 1:num_histories
    for mode_idx = 1:num_modes
        i = 2:n/2;
        nodal_matrix(i, mode_idx) = interp1(eigenfreqs, squeeze(modal_matrix(history_idx, mode_idx, :)), freqs(i));
        nodal_matrix(n/2+1, mode_idx) = real(interp1(eigenfreqs, squeeze(modal_matrix(history_idx, mode_idx, :)), freqs(n/2+1)));
        nodal_matrix(n/2+2:end, mode_idx) = conj(nodal_matrix(n/2:-1:2, mode_idx));
    end

    % Combine nodal contributions for all modes
    time_domain_signal(:, history_idx) = sum(nodal_matrix .* modal_process, 2);
end

%% 5. Perform inverse FFT to get time-domain samples
samples = zeros(n, num_histories);
for history_idx = 1:num_histories
    samples(:, history_idx) = ifft(time_domain_signal(:, history_idx));
end

end

function supports_vector = isVectorInputSupported(func_handle)
% Check if the function handle supports vector input
try
    func_handle([1, 2, 3, 4]); % Test function with a vector
    supports_vector = true;
catch
    supports_vector = false;
end
end
