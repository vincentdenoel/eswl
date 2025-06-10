function V_normalized = NormalizeMode(V,normMode)
% V_normalized = NormalizeMode(V)
% this function normalizes the matrix of modes V to a maximum real unitary
% for each mode
% V. Denoël, Nov. 2008

N = size(V,2);

if nargin==1; normMode =1; end

for s=1:N
    vec = V(:,s);
    veca = vec.*conj(vec);
    amax=veca(1);imax=1; for i=2:length(veca); if veca(i)>amax;amax=veca(i);imax=i;end;end
    
    if normMode == 1
        vec = vec / vec(imax);
    else
        vec = vec / norm(vec);
    end
    V_normalized(:,s) = vec;
end