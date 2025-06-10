function results = determineESWL(time,dp,nodalInfo,modalInfo)
% n  = number of measured time steps
% nb = number of pressure taps
% nl = number of loading points given by the design office
% nr = number of responses
% nm = number of modes for the dynamic analysis
%
% coord_taps : [nb x 3] (x,y,z)
% dp         : [n x nb] pressure time history
% nodes      : [nl x 3] (x,y,z)
% nodalInfo is a structure containing :
% .z         : [nl x nr] values of responses for unit loads at loading pts
% modalInfo is a structure containing :
% .z         : [nm x nr] values of responses for unit modal amplitude
% .freq      : [nm x 1] natural frequencies
% .mgen, .kgen : [nm x 1] modal masses and stiffnesses
% .damp      : [nm x 1] modal damping ratio (vs critical)
% .mode      : [nl x 3 x nm] mode shapes; values at each node (x,y,z)
% .mil       : [nl x 3 x nm] modal intertial loads; values at each node (x,y,z)

n  = size(dp,1); % number of time steps
nl = size(nodalInfo.z,1);
nr = size(nodalInfo.z,2);
nm = modalInfo.nm;

%n = round(n/15)  % Trim the signal to run faster. Do this for debug only.

% Time history of the pressure at nodes and QS responses
dp_nodes = dp'; % Forces at nodes

% Mean pressure field and fluctuations
mean_dp  = mean(dp_nodes,2);                  % mean net pressure at nodes
mean_F   = mean_dp;                           % mean force at nodes
fluct_dp = dp_nodes - repmat(mean_dp,1,n);    % fluct net pressure at nodes
fluct_F  = fluct_dp;                          % fluct force at nodes

% Covariance of wind forces at the nodes
cov_F = cov(fluct_dp');                    % covariance of forces at nodes
var_F = diag(cov_F);                       %  variances of forces at nodes
cov_F_coh = sqrt(var_F*var_F');

%% Background responses (nodal basis)
zn = nodalInfo.z;   % nl x nr : static response for unit loads at nodes
% Average responses 
mean_z = zn' * mean_F;
zb     = zn' * fluct_F;  % background component of the responses

% Fluctuations of responses and their covariances
fluct_z   = zn' * fluct_F;         % nr x n   Background responses over time
cov_z     = zn' * cov_F * zn ;     % Cov. of backg. Same result as cov(fluct_z');
cov_z_coh = zn' * cov_F_coh * zn ; % Cov. of backg. as if forces were fully coherent




%% Resonant responses (modal basis)
zm = modalInfo.z;   % nm x nr : responses for unit amplitudes in modes

F_star = zeros(nm,n);
for mode = 1:nm
    for node=1:nl
        F_star(mode,:) = F_star(mode,:) + fluct_F(node,:) * modalInfo.mode(node,2,mode);
    end
end

dt = time(2)-time(1);

for mode = 1:nm
    f = modalInfo.freq(mode);
    m = modalInfo.mgen(mode);
    xi= modalInfo.damp(mode);
    p = F_star(mode,:);
    [dep, ~, ~, dep_st] = Newmark(m,xi,f, p, dt,n-1,0.25,0.5);

    qb(mode,:) = dep_st;   % nm x n : Background response in each mode
    q(mode,:)  = dep;      % nm x n : Dynamic response in each mode

    qr(mode,:) = dep-dep_st;  % nm x n : Resonant part of the modal amplitude
    zr_mode(:,:,mode) = qr(mode,:)' * zm(mode,:);    % n x nr x nm
    std_zr_mode(:,mode) = std(zr_mode(:,:,mode)); % nr x nm  = std of zr_mode
end

zr = sum(zr_mode,3)'; % nr x n : Total resonant componant of the response

z  = zr + zb;  % n x nr = Sum of background and resonant components

% Determination of the standard deviations
std_q  = std(qr,[],2);
std_z  = std(z, [],2); % nr x 1
std_zb = std(zb,[],2); % nr x 1
std_zr = std(zr,[],2); % nr x 1


%% Do a basic extreme value statistics
T = 600;
nw = floor(T/dt);
nB = floor(n/nw);
nB = max([1 nB]);

for block=1:nB
    zminBlock(:,block)=min(z(:,1+(block-1)*nw:block*nw),[],2); % nB x nr
    zmaxBlock(:,block)=max(z(:,1+(block-1)*nw:block*nw),[],2); % nB x nr
end
zMin_mean = mean(zminBlock,2); % 1 x nr
zMax_mean = mean(zmaxBlock,2); % 1 x nr

gmin = zMin_mean./std_z; % 1 x nr
gmax = zMax_mean./std_z; % 1 x nr
g = (gmax-gmin)/2;       % we keep it symmetric in this illustration (*)
z_max = mean_z + g.*std_z;  % maximum response (upper envelope)
z_min = mean_z - g.*std_z;  % minimum response (lower envelope)

%(*) you can read out more on reconstruction of non-symmetric non-Gaussian
%    envelopes from my ORBI repository.

%% Computation of ESWL

eswl_b     = zeros(nl,nr);
eswl_b_coh = zeros(nl,nr);
eswl_r     = zeros(nl,nr,3);
ESWL_b     = zeros(nl,nr,3);
ESWL_r     = zeros(nl,nr,3);

wght_b = zeros( 1,nr);
wght_r = zeros(nm,nr);

z_rec = zeros(nr,nr);

% (*) see note at bottom of this routine for information about z_rec and ESWLs

for r = 1:nr

    % magnitude of force (normal to surface)
    eswl_b(:,r)     = g(r) * nodalInfo.z(:,r)' * cov_F     / std_zb(r);
    eswl_b_coh(:,r) = g(r) * nodalInfo.z(:,r)' * cov_F_coh / std_zb(r);

    wght_b(r) = std_zb(r) / std_z(r);
    
    % Components of the force in each direction (xyz)
    ESWL_b(:,r) = eswl_b(:,r) * wght_b(r);  % weighting
    
    % reconstructed responses (without average responses, in mean_z, to be added)
    z_rec(:,r) = (nodalInfo.z)' * eswl_b(:,r) * wght_b(r); %  background + ...
    
    for mode = 1:nm

        wght_r(mode,r) = std_zr_mode(r,mode)/std_z(r);

        % this is to make sure dynamic contributions to z_rec(r,r) are all positive
        wght_r(mode,r) = wght_r(mode,r) * sign(modalInfo.z(mode,r));
        
        eswl_r(:,r,1) = g(r) * std_q(mode) * modalInfo.mil(:,1,mode); % along x
        eswl_r(:,r,2) = g(r) * std_q(mode) * modalInfo.mil(:,2,mode); % along y
        eswl_r(:,r,3) = g(r) * std_q(mode) * modalInfo.mil(:,3,mode); % along z

        ESWL_r(:,r,:) = ESWL_r(:,r,:) + eswl_r(:,r,:) * wght_r(mode,r); % weighting

        z_rec(:,r) = z_rec(:,r) + modalInfo.z(mode,:)' * g(r) * std_q(mode) * wght_r(mode,r); % ... + resonant
    end

end

ESWL = ESWL_b + ESWL_r;

ESWL_m = mean_F; % mean forces


%% Information about z_rec and ESWLs
%
% z_rec(:,r) contains the structural responses obtained under the ESWL
% associated with response r (background and resonant fluctuation only)
%
% z_rec(r,r) is not exactly equal to g(r)*std_z(r) because it is assumed
% that resonant components are fully uncorrelated, which might not be
% totally right in case of close natural frequencies.
%
% The responses z_rec(:,r) are intended to be the most probable values
% of structural responses when response r reaches its max/min value.
%
% The responses z_rec(:,r) can be obtained under the following wind
% loading : ESWL(:,r,:) = ESWL_b(:,r,:) + ESWL_r(:,r,:)
%
% The min/max structural responses obtained under ESWL associated with
% response r are mean_z +/- z_rec(:,r). They correspond to the static loads
% ESWL_m(:,[xyz]) +/- ESWL(:,r,[xyz])

%% Prepare output



results.mean_z = mean_z;
results.std_z  = std_z;
results.std_zb = std_zb;
results.std_zr = std_zr;

results.wght_b = wght_b;
results.wght_r = wght_r;

results.z     = z;

results.z_min = z_min;
results.z_max = z_max;
results.z_rec = z_rec;

results.ESWL_m = ESWL_m;  % [nl x 3] : average wind forces at nodes in each
                          % direction (x,y,z) --> use these loads to
                          % recover mean_z
results.eswl_b = eswl_b;  % [nl x nr] : magnitude of ESWL for background part at each
                          % node and for each response
results.ESWL_b = ESWL_b;  % [nl x nr x 3] : 3 components of ESWL for background part
results.ESWL_r = ESWL_r;  % [nl x nr x 3] : 3 components of ESWL for resonant part
results.ESWL   = ESWL  ;  % sum of ESWL_b and ESWL_r


