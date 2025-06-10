function [results, ESWLresults] =  ...
    do_buffeting_analysis(windModel, structuralModel, aeroSec, simres_file)
% Do buffeting analysis of a FE model. In this simplified instance the
% geometry of the structure is a line and the height of all points is cst.
% This script sequentially :
% - run spectral analysis in frequency domain
% - run spectral analysis in time domain (for illustration purposes)
% - determines ESWLs for all 85 displacements and all 85 considered bending
%   moments 
%
% Input :
% - structuralModel  : mode shapes, daming, modal masses
% - aeroSec          : static and unsteady aerodynamic coeffs (CL,CM,Scanlan)
% - windModel        : turbulence (used for buffeting only)
% - simres_file      : string containing the file name where to save results

results = [];

%% interpret input

% Structural properties
M = structuralModel.M;
K = structuralModel.K;
C = structuralModel.C;

xi_s = structuralModel.xi_s;

XNOD   = structuralModel.XNOD;
NDOF   = structuralModel.NDOF;  % with only 2-dofs per node
nModes = structuralModel.nModes;
NElem  = length(structuralModel.ELEMNOA);

%% Compute mode shapes

[V, D]= eig(K,M);

V = V(:,1:nModes);
V = NormalizeMode(V);

omega = sqrt(diag(D));
omega = omega(1:nModes);
fnat  = omega/2/pi;

Mgen = diag(V'*M*V);
Kgen = diag(V'*K*V);
Cgen = 2*sqrt(Mgen(:).*Kgen(:)).*xi_s(:); % use damping ratio in modal basis

A_disp = V; % Influence matrix for displacements (=mode shapes)

for imod=1:nModes
    [~, ~, Bendi]=InternalForce(NDOF,V(:,imod),structuralModel.Ke,...
        structuralModel.corres,NElem,structuralModel.ELEMDOF);
    ABendi(:,imod) = [Bendi(:,1); Bendi(end,2)];
end

for idof=1:length(K)
    x_tmp = zeros(length(K),1);
    x_tmp(idof) = 1;
    [~, ~, Bendi]=InternalForce(NDOF,x_tmp,structuralModel.Ke,...
        structuralModel.corres,NElem,structuralModel.ELEMDOF);
    A_Bendi(:,idof) = [Bendi(:,1); Bendi(end,2)];
end

%% Spectral analysis in frequency domain

fprintf(' * Running spectral analysis in frequency domain...')

tmp = load('f_stoch'); % Read frequencies for which the spectral analysis is carried on
f_stoch = tmp.f_stoch; clear tmp

v = V(structuralModel.loadedDOFs,:); % keep loaded dofs only
kn = 1:size(v,1); % list of loaded nodes (1:85)

Sq = nan(nModes,nModes,length(f_stoch));
SqB= nan(nModes,nModes,length(f_stoch));
Sx = nan(NDOF,NDOF,length(f_stoch));
Sxf= nan(NDOF,length(kn),length(f_stoch));
Sf = nan(length(kn),length(kn),length(f_stoch));

for ifr = 1:length(f_stoch)

    w = 2*pi*f_stoch(ifr);

    % Create PSD matrix of wind loads
    x = XNOD(:);
    Sf_drag = Sf_model2(w, x, windModel, aeroSec);

    % PSD matrix of modal drag forces
    Sfgen = v'*Sf_drag(kn,kn)*v;

    H = FCTTransfertMDDL(diag(Mgen),diag(Kgen),diag(Cgen),w/2/pi);

    Sq(:,:,ifr) = H*Sfgen*H';
    SqB(:,:,ifr) = 1./diag(Kgen)*Sfgen*diag(1./Kgen);

    % PSD matrix of nodal drag forces
    H__ = FCTTransfertMDDL(M,K,C,w/2/pi);
    H = H__(:,structuralModel.loadedDOFs);

    Sx(:,:,ifr) = H*Sf_drag(kn,kn)*H';

    Sxf(:,:,ifr) = H*Sf_drag(kn,kn);

    Sf(:,:,ifr) = Sf_drag;
end

% Modal basis results
covq  = real(trapz(f_stoch,Sq,3));
covqB = real(trapz(f_stoch,SqB,3));
covx = A_disp*covq*A_disp';
stdx = sqrt(diag(covx));
covbending = ABendi*covq*ABendi';
stdbending = sqrt(diag(covbending));

% Nodal basis results
cov_f = real(trapz(f_stoch,Sf,3));
cov_x  = real(trapz(f_stoch,Sx,3));
cov_xf = real(trapz(f_stoch,Sxf,3));

std_x = sqrt(diag(cov_x));
cov_bending = A_Bendi*cov_x*A_Bendi';
std_bending = sqrt(diag(cov_bending));

results.f_stoch = f_stoch;
results.Sx = Sx;
results.Sxf = Sxf;
results.Sf = Sf;
results.Sq = Sq;

results.modal.covq = covq;
results.modal.cov_x = covx;
results.modal.std_x = stdx;
results.modal.cov_M = covbending;
results.modal.std_M = stdbending;

results.nodal.cov_x = cov_x;
results.nodal.cov_f = cov_f;
results.nodal.cov_xf = cov_xf;
results.nodal.std_x = std_x;
results.nodal.cov_M = cov_bending;
results.nodal.std_M = std_bending;

fprintf('Done.\n')

%% Spectral analysis in time domain, for illustration purposes
fprintf(' * Running spectral analysis in time domain...')
redo_generation = 0;

% File where to save/read the drag forces
dragforces_file = fullfile(cd, 'data', 'drag_forces.mat');

% If file does not exist create it.
if ~exist(dragforces_file, 'file')
    fprintf('Drag forces not found.\n')
    redo_generation = 1;
end


if ~redo_generation
    fprintf('Wind loads reloaded from previous run.\n')
    load drag_forces time drag_forces

else
    fprintf('Generating spectrum compatible wind loads...')
    dt = 0.04;
    S_hdl = @(w) Sf_model2(w, x, windModel, aeroSec)/4/pi;
    [time, drag_forces] = gen_samples(S_hdl, (0:2^16-1)*dt, [0 2.^(-9:0.05:4)]);
    save (fullfile('data','drag_forces.mat'), 'time', 'drag_forces')
    fprintf('Done.\n')
end

nl = size(drag_forces,2); % number of loading points

% unit loads at each 'pressure tap'
for i=1:nl
    p=zeros(length(K),1);
    p(2*i-1) = 1;
    x = K\p;
    [Axial, Shear, Bendi]=InternalForce(NDOF,x,structuralModel.Ke,...
        structuralModel.corres,NElem,structuralModel.ELEMDOF);
    nodalInfo.z_x(i,:) = x;
    nodalInfo.z_m(i,:) = [Bendi(:,1); 0];
end
nodalInfo.z   = [nodalInfo.z_x, nodalInfo.z_m];

mil = K*V;
modalInfo.z_x = V(:,:)';
modalInfo.z_m = ABendi';
modalInfo.z   = [modalInfo.z_x, modalInfo.z_m];
modalInfo.freq = fnat;
modalInfo.mgen = Mgen;
modalInfo.kgen = Kgen;
modalInfo.nm   = nModes;
modalInfo.damp = xi_s';
modalInfo.mode(:,1,:) = 0*V(1:2:end,:);
modalInfo.mode(:,2,:) = V(1:2:end,:);
modalInfo.mode(:,3,:) = 0;
modalInfo.mil(:,1,:) = 0*mil(1:2:end,:);
modalInfo.mil(:,2,:) = mil(1:2:end,:);
modalInfo.mil(:,3,:) = 0;

ESWLresults = determineESWL(time,drag_forces,nodalInfo,modalInfo);

save(simres_file, 'results', 'ESWLresults', 'modalInfo', 'nodalInfo', ...
    'time', 'drag_forces', 'structuralModel', 'windModel', 'aeroSec')

fprintf ('Done.')
