function [x_rec, m_rec, scale] = plot_diagrams(ESWLresults,structuralModel,i_x,i_m,i_z,clr,eswl, fit_env_scale, doplot)

% i_x = indices in z with displacements and rotations of nodes
% i_m = indices in z with bending moments
% i_z = list of responses

if nargin<8, fit_env_scale = 0; end
if nargin<9, doplot = 1; end

XNOD = structuralModel.XNOD;

upper = ESWLresults.z_max;
lower = ESWLresults.z_min;

if nargin>=7
    % Recompute the displacement and bending moment for a given eswl (and scale
    % it)

    NDOF = structuralModel.NDOF;  % with only 2-dofs per node
    K = structuralModel.K;

    %eswl = ESWLresults.ESWL_b(:,:,2);
    f = zeros(size(eswl,1)*2,length(i_z));
    f(1:2:end,:)=eswl(:,i_z);
    xx = K\f;
    
    [~, ~, Bendi]=InternalForce(NDOF,xx,structuralModel.Ke,...
        structuralModel.corres,length(structuralModel.ELEMNOA),structuralModel.ELEMDOF);
    

    if fit_env_scale
        scale_x = min(abs( upper(i_x)./xx(1:2:end) )); % make sure there is no overestimation
        scale_m = min(abs( upper(i_m(1:end-1))./Bendi(:,1) ));
        scale = min([scale_x scale_m]);
    else
        scale = upper(i_z) / xx(i_z); % match envelopes at point where response is recovered (work if i_z is a displacement) - adapt code otherwise
        %scale = 1;
    end
    xx =   scale  * xx;
    Bendi= scale * Bendi;
    eswl = scale * eswl;

    x_rec = xx(1:2:end,:);

    m_rec = [Bendi(:,1); Bendi(end,2)];
    
else
    eswl = ESWLresults.ESWL(:,:,2);
    x_rec= ESWLresults.z_rec(i_x,i_z);
    m_rec= ESWLresults.z_rec(i_m,i_z);
    scale = 1;
end

if ~doplot, return, end


if i_z<=170
    i__z = (i_z+1)/2; % this is a displacement
else
    i__z = i_z-170; % this is a bending moment
end

plotSym = 1;

subplot(311)
plot(XNOD,upper(i_x),'k:','LineWidth',1.2), hold on
plot(XNOD,lower(i_x),'k:','LineWidth',1.2)
plot(XNOD,x_rec,'color',clr)
if plotSym, plot(XNOD,-x_rec,'color',clr), end
plot(XNOD(i__z),x_rec(i__z),'ko','MarkerFAcecolor','k')
title ('Displacements')
xlim tight, box off

subplot(312)
plot(XNOD,upper(i_m),'k:','LineWidth',1.2), hold on
plot(XNOD,lower(i_m),'k:','LineWidth',1.2)
plot(XNOD,m_rec,'color',clr)
if plotSym, plot(XNOD,-m_rec,'color',clr), end
plot(XNOD(i__z),m_rec(i__z),'ko','MarkerFAcecolor','k')
title ('Bending moment')
xlim tight, box off


subplot(313)
plot(XNOD,eswl(:,i_z),'color',clr), xlim tight, box off,  hold on
plot(xlim, [0 0],'k','linewidth',0.3)
title ('Static wind loading')

set(gcf,'position',[ 218   433   782   364])
pause(0.3)