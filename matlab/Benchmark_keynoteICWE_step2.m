% MATLAB code to produce the figures in:
% V. Denoël, Quantifying Complexity in the Envelope Reconstruction Problem:
% Review, Comparison, and a Detailed Illustration, 
% Journal of Wind Engineering & Industrial Aerodynamics, 2024
%
% Please reference the above paper if you find this script useful.

%% Convergence plots for combinaisions of PSWL

close all


ifig = figure;
target_x = abs(ESWLresults.z_max(i_x));
target_m = abs(ESWLresults.z_max(i_m));

i_z=1;
clear x_rec
clear m_rec
clear reconst_rate

x_rec = [];
m_rec = [];

nb_comb = 14; % number of combinations

for i_comb=1:nb_comb

    figure(ifig)
    iclr=i_comb; if iclr>7;iclr=iclr-7; end
    clr = cmp(iclr,:);


    x_rec_trial = x_rec;
    m_rec_trial = m_rec;

    for i_trial = 1:20000
        
        trial_coefs(:,i_trial,i_comb) = randn(7,1);
        trial_eswl = U(:,1:7)*trial_coefs(:,i_trial,i_comb);
        [x_rec_trial(:,i_comb), m_rec_trial(:,i_comb)] = plot_diagrams(ESWLresults,structuralModel,i_x,i_m,i_z,clr,trial_eswl, 1, 0);

        x_env_trial = max(1.*[x_rec_trial -x_rec_trial],[],2);
        m_env_trial = max(1.*[m_rec_trial -m_rec_trial],[],2);

        x_env_tmp = x_env_trial;
        m_env_tmp = m_env_trial;
        ichk_x = find(x_env_trial>target_x);
        ichk_m = find(m_env_trial>target_m);
        x_env_tmp(ichk_x) = target_x(ichk_x);
        m_env_tmp(ichk_m) = target_m(ichk_m);
        reconst_rate_trial(i_trial,:) = [sum(x_env_tmp)/metric_x sum(m_env_tmp)/metric_m ];
    end
    % keep best test
    [~, i_t] = max(mean(reconst_rate_trial'));

    coefs(:,i_comb) = trial_coefs(:,i_t,i_comb);
    reconst_rate(i_comb,1:2) = reconst_rate_trial(i_t,1:2); 


    eswl = U(:,1:7)*coefs(:,i_comb);
    [x_rec(:,i_comb), m_rec(:,i_comb)] = plot_diagrams(ESWLresults,structuralModel,i_x,i_m,i_z,clr,eswl, 1);
    x_env = max(1.*[x_rec -x_rec],[],2);
    m_env = max(1.*[m_rec -m_rec],[],2);



    figure(999),
    plot(reconst_rate), grid, ylim ([0.3 1]), xlim([0 10])
    set(gcf,'position',[559   593   441   204])
    

    figure(100)
    sub1(i_comb)=subplot(14,2,2*i_comb-1); plot(x_env), hold on, plot(-x_env), xlim tight
    plot(abs(ESWLresults.z_max(i_x))),plot(-abs(ESWLresults.z_max(i_x)))
    sub2(i_comb)=subplot(14,2,2*i_comb  ); plot(m_env), hold on, plot(-m_env), xlim tight
    plot(abs(ESWLresults.z_max(i_m))),plot(-abs(ESWLresults.z_max(i_m)))    
    drawnow
end


figure(100),

ax1 = axis(subplot(14,2,27));
ax2 = axis(subplot(14,2,28));
for i_comb=1:14
    axes(sub1(i_comb)), axis(ax1), axis off, set(gca,'position',[0.01 0.92-(i_comb-1)*0.07 0.48  0.07])
    axes(sub2(i_comb)), axis(ax2), axis off, set(gca,'position',[0.51 0.92-(i_comb-1)*0.07 0.48  0.07])
end
set(gcf,'Position',[1050          81         378         716])

%save reconstruction


