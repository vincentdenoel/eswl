% MATLAB code to produce the figures in:
% V. Denoël, Quantifying Complexity in the Envelope Reconstruction Problem:
% Review, Comparison, and a Detailed Illustration, 
% Journal of Wind Engineering & Industrial Aerodynamics, 2024
%
% Please reference the above paper if you find this script useful.

clc;
clear;
close all;

% Add custom tools for structural and stochastic dynamics, and paths to data
addpath('tools', 'data');

% Load finite element model
% Expected variables:
%   - Structural matrices: K, M, C
%   - Element matrices: Ke
%   - Node and element data: NNode, NElem, XNOD, ELEMNOA, ELEMNOB
%   - DOF data: ELEMDOF, CORRES, NDOF, iDOF_obs
load finite_element_model;

% Redo analysis: Set to 1 to recalculate and overwrite `SimResults.mat`
redo_analysis = 0;

% File where to save/read the data for the envelope reconstruction analysis
simres_file = fullfile(cd, 'data', 'SimResults.mat');

% If file does not exist create it.
if ~exist(simres_file, 'file')
    fprintf('Results from previous analysis not found.\n')
    redo_analysis = 1;
end

if redo_analysis
    fprintf('Generation of wind buffeting analysis data...\n')
    % Define wind model for buffeting analysis
    windModel.rho   = 1.22;      % Air density [kg/m³]
    windModel.U     = 34.66;     % Mean wind speed [m/s]
    windModel.sigu  = 4.56;      % Wind speed standard deviation
    windModel.Iu    = 4.56 / 34.66;  % Turbulence intensity
    windModel.Lux   = 50;        % Longitudinal turbulence scale [m]
    windModel.C     = 8;         % Additional aerodynamic parameter

    % Define aerodynamic section properties
    aeroSec.B  = 30;             % Deck width [m]
    aeroSec.CD = 0.4;            % Drag coefficient at zero AoA

    % Define structural model properties
    structuralModel.NDOF = NDOF;
    structuralModel.M = M;
    structuralModel.K = K;
    structuralModel.C = C;
    structuralModel.Ke = Ke;
    structuralModel.nModes = 7;
    structuralModel.xi_s = 0.003 * ones(1, 7); % Modal damping ratios

    structuralModel.XNOD = XNOD;
    structuralModel.ELEMNOA = ELEMNOA;
    structuralModel.ELEMNOB = ELEMNOB;
    structuralModel.ELEMDOF = ELEMDOF;

    structuralModel.iDOF_obs = iDOF_obs;
    structuralModel.corres = CORRES;
    structuralModel.loadedDOFs = find(mod(CORRES, 3)); % Identify loaded DOFs

    % Perform buffeting analysis and save results
    do_buffeting_analysis(windModel, structuralModel, aeroSec, simres_file);

end
