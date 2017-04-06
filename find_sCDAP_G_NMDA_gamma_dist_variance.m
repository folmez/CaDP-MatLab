function varargout = find_sCDAP_G_NMDA_gamma_dist_variance(varargin)
% This routine is to find the variance of the gamma distribution for a
% given spike timing difference so that the resulting calcium peak
% coefficient of variation is similar to what is predicted by the
% theoretical analysis (Figure 1D in Shouval 2004).

% Input arguments
spike_timing_diff = varargin{1};    % t_post - t_pre
nr_sims           = 250;
CV_rel_tol        = 0.01;
display_summary   = 1;

% Calculate theoretical CV as described in Shouval 2004
CV_theory         = 0.095 + 0.0045*max([0 spike_timing_diff]) - ...
    0.00067*min([0 spike_timing_diff]);

% Display summary header
if display_summary
    fprintf('\t SPIKE TIMING DIFFERENCE = %1.2f\n', spike_timing_diff);
    fprintf('Time\tG_NMDA gamma dist. var.\tCV sim\tCV theory\t');
    fprintf('CV rel error (tolerance = %1.2f)\n', CV_rel_tol);
    [~, tSIM] = display_sim_progress('initialize');
end

% Model
sigma_square_small   = 1e-8;
sigma_square_large   = 1e-4;
while 1
    % Calculate next variance to be tested
    sigma_square_next = sqrt(sigma_square_small*sigma_square_large);
    
    % Find CV corresponding to the next 
    CV_sim = find_CV_sim(spike_timing_diff, sigma_square_next, nr_sims);
    
    % Find absolute error and relative error of CV
    CV_abs_error = abs(CV_sim-CV_theory);
    CV_rel_error = CV_abs_error/CV_theory;
    
    % Display summary
    if display_summary
        fprintf('[%1.2fm]\t', toc(tSIM)/60);
        fprintf('\t%1.4e\t', sigma_square_next);
        fprintf('%1.4f\t', CV_sim);
        fprintf('%1.4f\t', CV_theory);
        fprintf('\t%1.4f\n', CV_rel_error);
    end
    
    if CV_rel_error < CV_rel_tol
        % break out of loop is relative CV tolerance is achieved
        break;
    elseif CV_sim > CV_theory
        % update large gamma distribution variance if sim CV is bigger
        sigma_square_large = sigma_square_next;
    elseif CV_sim < CV_theory
        % update small gamma distribution variance if sim CV is smaller
        sigma_square_small = sigma_square_next;
    end
    
end

% Output arguments
varargout{1} = sigma_square_next;

end


%%
function CV_sim = find_CV_sim(varargin)

% Input arguments
spike_timing_diff = varargin{1};    % t_post - t_pre
sigma_square      = varargin{2};
nr_sims           = varargin{3};

NMDAr_version     = 'stochastic';   
display_progress  = 0;
plot_results      = 0;
Z                 = 10;           % # on NMDA receptors
dt                = 1;              % in ms
% --------------------------------------------------------------------
i = 6;
while i<=length(varargin),
    switch varargin{i},
        case 'Z',                              Z = varargin{i+1};
        case 'dt',                            dt = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% --------------------------------------------------------------------
BPAP_delay   = 0;

t_pre_spike    = 0;                 % Pre-synaptic spike time
t_post_spike   = spike_timing_diff; % Postsynaptic spike time

closed_NMDAr_frac_before_spike = 1;

V_rest         = -65;               % Resting potential
EPSP_norm      = 0.6968;            % See Shouval 2002 supplementary
EPSP_param_set = 1;
t0             = min(t_pre_spike, t_post_spike);
tend           = 500;
nr_time_steps  = (tend-t0)/dt+1;

% Initialize
Ca_peaks = zeros(nr_sims, 1);

% Model
[tDISP, tSIM] = display_sim_progress('initialize');
for k = 1:nr_sims
    % Initialize
    V_post = V_rest * ones(nr_time_steps, 1);   % Postsyn. membrane pot.
    NMDAr_frac = zeros(nr_time_steps, 1);       % Open NMDAr frac
    NMDAr_cal_cur = zeros(nr_time_steps, 1);    % NMDAr calcium current
    Ca = zeros(nr_time_steps, 1);               % Calcium level

    % No need to consider the case when G_NMDA is 0 because later CV is
    % calculated from nonzero Ca peaks. So here we can force G_NMDA to be
    % nonzero. This should cut simulation time in half without losing
    % accuracy
    while 1
        G_NMDA = generate_random_maximal_chord_conductance(...
            spike_timing_diff, Z, 'sigma_square', sigma_square);
        if G_NMDA ~= 0
            break;
        end
    end

    % Model
    for i = 2:nr_time_steps
        t_next = t0+i*dt;
        % Calculate potential at the dentrite
        V_post(i) = V_rest + EPSP(t_next, t_pre_spike, EPSP_norm, ...
            EPSP_param_set) + ...
            BPAP(t_next, t_post_spike, 'version', 'stochastic', ...
            'BPAP_delay', BPAP_delay);
        % Calculate NMDA calcium current
        [NMDAr_cal_cur(i), NMDAr_frac(i)] = ...
            NMDAr_calcium_current(t_next, t_pre_spike, V_post(i), ...
            closed_NMDAr_frac_before_spike, ...
            'version', NMDAr_version, ...
            'G_NMDA_for_each_spike', G_NMDA);
        % Update Calcium level
        Ca(i) = update_Ca(Ca(i-1), NMDAr_cal_cur(i), dt);
    end
    
    % Record Calcium peak
    Ca_peaks(k) = max(Ca);
    
    % Display simulation progress
    tDISP = display_sim_progress(k, nr_sims, tDISP, tSIM, display_progress);
end

% Calculate simulation CV from nonzero Ca peaks
CV_sim    = std(nonzeros(Ca_peaks))/mean(nonzeros(Ca_peaks));

% Plot results
if plot_results
figure,
    if nr_sims == 1
        tt = linspace(t0, tend, nr_time_steps)';
        plot(tt, Ca);
        h_leg = legend('Calcium level', 'Location', 'Best');
        set(h_leg, 'FontSize', 15);
        xlabel('Time (ms)', 'FontSize', 15);
    elseif nr_sims > 1
        CV_abs_error  = CV_sim - CV_theory;
        histogram(Ca_peaks, 'Normalization', 'probability');
        title(['\Deltat = ' num2str(spike_timing_diff) ...
            ', nonzero Ca peaks CV = ', ...
            num2str(CV_sim, '%1.2f'), ...
            ', sim-theory =  ', num2str(CV_abs_error, '%1.2f')]);
    end
end

end
