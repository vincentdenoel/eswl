% MATLAB code to produce the figures in:
% V. Denoël, Quantifying Complexity in the Envelope Reconstruction Problem:
% Review, Comparison, and a Detailed Illustration,
% Journal of Wind Engineering & Industrial Aerodynamics, 2024
%
% Please reference the above paper if you find this script useful.

clear, close all

% add some of my tools for structural and stochastic dynamics, add path to
% additional data
addpath ('tools', 'data')

cmp = get(groot, 'DefaultAxesColorOrder');

% Simulation results are generated in step 0. Recover them now.
simres_file = fullfile(cd, 'data', 'SimResults.mat');
if ~exist(simres_file, 'file')
    run('Benchmark_keynoteICWE_step0.m')
else
    load (simres_file)
end

% Recover and store variables (from structure 'results')
fieldNames = fieldnames(results);

for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    fieldValue = results.(fieldName);
    eval([fieldName ' = fieldValue;']);
end

XNOD =     structuralModel.XNOD;
%iDOF_obs = structuralModel.iDOF_obs;
%% plot PSD of some responses

iDOF_obs = structuralModel.iDOF_obs([2 4 6]);

figure,
plotTF(time, ESWLresults.z(iDOF_obs, :), 2048*4)
xlim ([0 2])
title ('PSD of some displacements (middle of 2nd, 4th, 6th spans)')


% plot stochastic spectral analysis for comparison
figure
for i=iDOF_obs
    semilogy(f_stoch, squeeze(Sx(i, i, :)))
    xlim ([0 2]), hold on
end
title ('PSD of some displacements (middle of 2nd, 4th, 6th spans)')


%% Compute and plot CPT modes
i_x = 1:2:170;
i_m = 171:255;

[Vf, Df] = eigs(nodal.cov_f,7);
figure,
for i=1:7
    subplot(7, 1, i); plot(Vf(:, i)); axis tight;axis off;
    title (Df(i,i))
    hold on; plot(xlim, [0 0], 'k')
end

close all
for imode = [1 3 5]

    i_z = structuralModel.iDOF_obs(2);

    eswl(:, i_z) = Vf(:, imode);
    %eswl(:, i_z) = ones(size(Vf(:, 1)));
    clr = cmp(imode, :);
    plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, eswl)
end

%% LRC - displacements

iDOF_obs = structuralModel.iDOF_obs([2 4 6]);

f_LRC = zeros(length(structuralModel.corres), length(iDOF_obs));
x_reconstruct = zeros(length(structuralModel.corres), length(iDOF_obs));

% from time series
enveloppe = std(ESWLresults.z(1:2:170, :)');
for i=1:length(iDOF_obs)
    cov_xf = (ESWLresults.z(iDOF_obs(i), :) * drag_forces(:, :))/65536;
    std_x =   std(ESWLresults.z(iDOF_obs(i), :));

    f_LRC(structuralModel.loadedDOFs, i) = cov_xf / std_x;
    x_reconstruct(:, i) = structuralModel.K \ f_LRC(:, i);

    scale = x_reconstruct(iDOF_obs(i), i) / std_x;

    f_LRC(structuralModel.loadedDOFs, i) = 1/scale * f_LRC(structuralModel.loadedDOFs, i);
    x_reconstruct(:, i) = 1/scale * x_reconstruct(:, i);
end

figure
for i=1:length(iDOF_obs)
    i_z = iDOF_obs(i);
    eswl(:, i_z) = f_LRC(1:2:end, i);

    clr = cmp(i, :);
    plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, eswl)
end

% from spectral analysis
enveloppe = nodal.std_x(1:2:end);
for i=1:length(iDOF_obs)

    cov_xf = nodal.cov_xf(iDOF_obs(i), :);
    std_x  = nodal.std_x (iDOF_obs(i));

    f_LRC(structuralModel.loadedDOFs, i) = cov_xf / std_x;
    x_reconstruct(:, i) = structuralModel.K \ f_LRC(:, i);

    scale = x_reconstruct(iDOF_obs(i), i) / std_x;


    f_LRC(structuralModel.loadedDOFs, i) = 1/scale * f_LRC(structuralModel.loadedDOFs, i);
    x_reconstruct(:, i) = 1/scale * x_reconstruct(:, i);
end

figure
for i=1:length(iDOF_obs)

    i_z = iDOF_obs(i);
    eswl(:, i_z) = f_LRC(1:2:end, i);

    clr = cmp(i, :);
    plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, eswl)
end



%% ESWL -- Plot a few of them
cmp = get(groot, 'DefaultAxesColorOrder');
close all

i_x = 1:2:170;
i_m = 171:255;
i_z = structuralModel.iDOF_obs([2 4 6]);

clr = cmp(1, :);
plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr)


figure
i_z = structuralModel.iDOF_obs(2);
i_z = structuralModel.iDOF_obs(1);
%i_z = 183 + [0:12:60];
clr = cmp(1, :);
plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr)
%%
close all
%i_z = structuralModel.iDOF_obs(1)-10;
i_z = structuralModel.iDOF_obs(5)+10;

eswl = ESWLresults.ESWL(:, :, 2);
clr = cmp(1, :);
plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, eswl)

eswl = ESWLresults.ESWL_b(:, :, 2);
clr = cmp(1, :);
plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, eswl)

eswl = ESWLresults.ESWL_r(:, :, 2);
clr = cmp(3, :);
plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, eswl)


%% PSWL
eswl = ESWLresults.ESWL(:, :, 2);
eswl(:, end)=2*eswl(:, end-1)-eswl(:, end-2); % compensate for NaN's in the last column (moment on support)

figure
subplot(121), contourf(eswl(:, i_x)'), axis tight, box off
xlabel('Abscissa x [#]'), ylabel ('ESWL - Displacement [#]')
subplot(122), contourf(eswl(:, i_m)'), axis tight, box off
xlabel('Abscissa x [#]'), ylabel ('ESWL - Moment [#]')

[U, S, V] = svd(eswl);
figure
bar(diag(S))
title ('Eigenvalue - Aerodynamic-Structural Complexity = 7')
ASC = 7;

figure
for i=1:ASC
    subplot(ASC, 1, i)
    plot(U(:, i)), axis tight, hold on
    plot(xlim, [0 0], 'k'), axis off
    title (sprintf('Principal static wind load #%d', i))
end
%% Plot response diagram under PSWLs
% They are scaled, there's no overshooting
close all
i_z = 1;

for i_pswl=1:3
    iclr=i_pswl; if iclr>7;iclr=iclr-7; end
    clr = cmp(iclr, :);
    [x_rec(:, i_pswl), m_rec(:, i_pswl)] = plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, U(:, i_pswl), 1);
end
x_env = max([x_rec -x_rec], [], 2);
m_env = max([m_rec -m_rec], [], 2);

%% Here we define a metric to quantify the envelope reconstruction rate
metric_x = sum(abs(ESWLresults.z_max(i_x)));
metric_m = sum(abs(ESWLresults.z_max(i_m)));

%% Convergence plots PSWL
fig = figure;
target_x = abs(ESWLresults.z_max(i_x));
target_m = abs(ESWLresults.z_max(i_m));
i_z = 1;
clear x_rec
clear m_rec
clear reconst_rate
for i_pswl=1:ASC
    figure(fig)
    iclr=i_pswl; if iclr>7;iclr=iclr-7; end
    clr = cmp(iclr, :);
    [x_rec(:, i_pswl), m_rec(:, i_pswl)] = plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, U(:, i_pswl), 1);
    x_env = max(1.*[x_rec -x_rec], [], 2);
    m_env = max(1.*[m_rec -m_rec], [], 2);

    x_env_tmp = x_env;
    m_env_tmp = m_env;
    ichk_x = find(x_env>target_x);
    ichk_m = find(m_env>target_m);
    x_env_tmp(ichk_x) = target_x(ichk_x);
    m_env_tmp(ichk_m) = target_m(ichk_m);
    reconst_rate(i_pswl, :) = [sum(x_env_tmp)/metric_x sum(m_env_tmp)/metric_m ];

    pause(0.5)
end
figure(999), plot(reconst_rate), grid, ylim ([0.3 1]), xlim([0 14]), hold on
set(gcf, 'position', [559   593   441   204])

%% Convergence plots Naive approach
figure
target_x = abs(ESWLresults.z_max(i_x));
target_m = abs(ESWLresults.z_max(i_m));

clear x_rec
clear m_rec
clear reconst_rate

eswl = ESWLresults.ESWL(:, :, 2);

for icas=1:14
    iclr=icas; if iclr>7;iclr=iclr-7; end

    if icas<=7
        i_z = structuralModel.iDOF_obs(icas);
    else
        i_z = structuralModel.iDOF_obs(icas-7)+12;
    end

    clr = cmp(iclr, :);
    [x_rec(:, icas), m_rec(:, icas)] = plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, eswl, 1);
    x_env = max(1.*[x_rec -x_rec], [], 2);
    m_env = max(1.*[m_rec -m_rec], [], 2);

    x_env_tmp = x_env;
    m_env_tmp = m_env;
    ichk_x = find(x_env>target_x);
    ichk_m = find(m_env>target_m);
    x_env_tmp(ichk_x) = target_x(ichk_x);
    m_env_tmp(ichk_m) = target_m(ichk_m);
    reconst_rate(icas, :) = [sum(x_env_tmp)/metric_x sum(m_env_tmp)/metric_m ];

    pause(0.5)
end

figure(999), plot(reconst_rate), grid, ylim ([0.3 1]), xlim([0 10]), grid on
set(gcf, 'position', [559   593   441   204])

%% Convergence plots Fastest Descent
%close all
fig = figure;

target_x = abs(ESWLresults.z_max(i_x));
target_m = abs(ESWLresults.z_max(i_m));

clear x_rec
clear m_rec
clear reconst_rate

eswl = ESWLresults.ESWL(:, :, 2);

i_z = structuralModel.iDOF_obs(1);

for icas=1:14
    figure(fig)
    iclr=icas; if iclr>7;iclr=iclr-7; end


    clr = cmp(iclr, :);
    [x_rec(:, icas), m_rec(:, icas)] = plot_diagrams(ESWLresults, structuralModel, i_x, i_m, i_z, clr, eswl, 1);
    x_env = max(1.*[x_rec -x_rec], [], 2);
    m_env = max(1.*[m_rec -m_rec], [], 2);

    x_env_tmp = x_env;
    m_env_tmp = m_env;
    ichk_x = find(x_env>target_x);
    ichk_m = find(m_env>target_m);
    x_env_tmp(ichk_x) = target_x(ichk_x);
    m_env_tmp(ichk_m) = target_m(ichk_m);
    reconst_rate(icas, :) = [sum(x_env_tmp)/metric_x sum(m_env_tmp)/metric_m ];

    [~, i_xx]= max(target_x - x_env_tmp);
    [~, i_mm]= max(target_m - m_env_tmp);
    if max(target_x - x_env_tmp) / metric_x > max(target_m - m_env_tmp) / metric_m
        i_z = 2*i_xx-1;
    else
        i_z = 170+i_mm;
    end


    figure(101)
    sub1(icas)=subplot(14, 2, 2*icas-1); plot(x_env), hold on, plot(-x_env), xlim tight
    plot(abs(ESWLresults.z_max(i_x))), plot(-abs(ESWLresults.z_max(i_x)))
    sub2(icas)=subplot(14, 2, 2*icas  ); plot(m_env), hold on, plot(-m_env), xlim tight
    plot(abs(ESWLresults.z_max(i_m))), plot(-abs(ESWLresults.z_max(i_m)))
    drawnow

end
figure(999), plot(reconst_rate), grid, ylim ([0.3 1]), xlim([0 14]), hold on
set(gcf, 'position', [559   593   441   204])


figure(101)
ax1 = axis(subplot(14, 2, 27));
ax2 = axis(subplot(14, 2, 28));
for icas=1:14
    axes(sub1(icas)), axis(ax1), axis off, set(gca, 'position', [0.01 0.92-(icas-1)*0.07 0.48  0.07])
    axes(sub2(icas)), axis(ax2), axis off, set(gca, 'position', [0.51 0.92-(icas-1)*0.07 0.48  0.07])
end