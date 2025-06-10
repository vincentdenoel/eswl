function Sf_drag = Sf_model2(omega, x, windModel, aeroSec)
%
% Returns the x-PSD matrix of drag loads
%
% windModel.U   : mean wind velocity
% windModel.Iu  : turbulence intensity u(t)
% windModel.Lux : turbulence lengthscale u,x
%

% Compute xPSD of turbulence components
nx = length(x);
n  = length(omega);

Sf_drag = zeros(nx,nx,n);

rho= windModel.rho; % kg/m3
U  = windModel.U;   % m/s
Iu = windModel.Iu;
L  = windModel.Lux;
sig_u = Iu*U;

B  = aeroSec.B;     % m

q = 1/2*rho*U^2;     % N/m2
C =  windModel.C;    % -

dx = abs(ones(nx,1)*x' - x*ones(1,nx)); % matrix of distances between nodes

% assume that the wind PSD is the same everywhere
f  = omega/2/pi;
y  = f*L/U; % Monin frequency
Su = 4*L/U*sig_u^2 ./ (1+70.7*y.^2).^(5/6);  % m^2/s

% Aerodynamic properties, should be taken from aeroSec
coef_D_u = q*B * (2*aeroSec.CD/U);    % Ns/m^2

% Consider exponential coherence
for i=1:n
    coh = exp (-C*omega(i)*dx/2/pi/U);
    Sf_drag(:,:,i)    = coef_D_u^2 * Su(i) * coh;   % N^2s/m^2
end

% Multiply by the length of each element
dx2 = (x(2:end)-x(1:end-1))/2;
DX  = ([dx2;0]+[0;dx2]);
DX  = DX*DX';
for i=1:n
    Sf_drag(:,:,i)    = Sf_drag(:,:,i)  .* DX;   % N^2s
end
