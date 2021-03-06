function varargout = CaDP(varargin)
sim_protocol = varargin{1};

%% SINGLE PRE, SINGLE PRE-POST, SINGLE POST-PRE
if strcmp(sim_protocol, 'Single')
    % Reference: Figure 2 in Shouval 2002
    % Figure 2A12: CaDP('Single', [-50, 250, -1e-5, inf], 'Narrow BPAP');
    % Figure 2A34: CaDP('Single', [-50, 250, -1e-5, -1e-5-10], 'Narrow BPAP');
    % Figure 2A56: CaDP('Single', [-50, 250, -1e-5, -1e-5+10], 'Narrow BPAP');
    % Figure 2B12: CaDP('Single', [-50, 250, -1e-5, inf], 'BPAP + ADP');
    % Figure 2B34: CaDP('Single', [-50, 250, -1e-5, -1e-5-10], 'BPAP + ADP');
    % Figure 2B56: CaDP('Single', [-50, 250, -1e-5, -1e-5+10], 'BPAP + ADP');
    
    % Input arguments
    t0           = varargin{2}(1);   % in ms
    tend         = varargin{2}(2);   % in ms
    t_pre_spike  = varargin{2}(3);   % in ms
    t_post_spike = varargin{2}(4);   % in ms
    BPAP_type    = varargin{3};
    NMDAr_I_f    = 0.50;             % NMDAr fast decay component 
    dt           = 0.1;              % in ms
    V_rest       = -65;              % in mV
    plot_NMDAr_calcium_current = 0;
    
    nr_time_steps = (tend-t0)/dt+1;
    
    % Initialization
    V_post = V_rest * ones(nr_time_steps, 1);   % Postsyn. membrane pot.
    NMDAr_cal_cur = zeros(nr_time_steps, 1);    % NMDAr calcium current
    NMDAr_frac = zeros(nr_time_steps, 1);       % Open NMDAr fraction
    Ca = zeros(nr_time_steps, 1);               % Calcium level
    closed_NMDAr_frac_before_spike = 1;        % Fraction of closed NMDAr
                                                % before the single spike
    
    % Model
    for i = 2:nr_time_steps
        t_next = t0+i*dt;
        % Calculate potential at the dentrite
        V_post(i) = V_rest + EPSP(t_next, t_pre_spike, 1, 1) + ...
            BPAP(t_next, t_post_spike, 'BPAP_type', BPAP_type);
        % Calculate NMDA calcium current
        [NMDAr_cal_cur(i), NMDAr_frac(i)] = NMDAr_calcium_current( ...
            t_next, t_pre_spike, V_post(i), ...
            closed_NMDAr_frac_before_spike, ...
            'NMDAr_I_f', NMDAr_I_f);
        % Update Calcium level
        Ca(i) = update_Ca(Ca(i-1), NMDAr_cal_cur(i), dt);
    end
    
    % Plot results
    tt = linspace(t0, tend, nr_time_steps)';
    figure,
    subplot(2,1,1)
    AX = plotyy(tt, NMDAr_frac, tt, V_post);
    h_leg = legend('NMDAr fraction', 'Post potential', ...
        'Location', 'Best');
    set(h_leg, 'FontSize', 15);
    set(gca, 'XTick', []);
    if t_post_spike == inf
        plot_title = 'Only pre';
    else
        plot_title = ['\Deltat = ' num2str(t_post_spike-t_pre_spike)];
    end
    title(plot_title, 'FontSize', 15);
    AX(1).YLim = [0    1];
    AX(2).YLim = [-70 35];
    
    subplot(2,1,2)
    plot(tt, Ca);
    ylim([0 0.8]);
    h_leg = legend('Calcium level', 'Location', 'Best');
    set(h_leg, 'FontSize', 15);
    xlabel('Time (ms)', 'FontSize', 15);
    hold on;
    plot([t0 tend], [0.35 0.35], 'k--');
    plot([t0 tend], [0.55 0.55], 'k--');
    
    if plot_NMDAr_calcium_current
        figure,
        plot(tt, NMDAr_cal_cur);
        h_leg = legend('NMDAr Ca current', 'Location', 'Best');
        set(h_leg, 'FontSize', 15);
        xlabel('Time (ms)', 'FontSize', 15);
    end
end

%% PAIRING PRESYNAPTIC STIMULATION WITH POSTSYNAPTIC VOLTAGE CLAMP
if strcmp(sim_protocol, 'PPSwPVC')
    % Reference: Figure 3A in Shouval 2002
    % Figure 3A: CaDP('PPSwPVC', [0 1e5], [1 100], [-65:1:-10], 0.5);
    
    % Input arguments
    t0             = varargin{2}(1);    % in ms
    tend           = varargin{2}(2);    % in ms
    pre_spike_freq = varargin{3}(1);    % Presynaptic freq (in Hz)
    nr_pre_spikes  = varargin{3}(2);    % # of presynaptic pulses
    V_post_c_vec   = varargin{4};       % Fixed post. potential (in mV)
    NMDAr_I_f     = varargin{5};       % NMDAr fast decay component 
                                        % relative magnitude
    
    dt             = 1;                 % in ms
    plot_only_w    = 1;
    w_init         = 0.25;              % Initial weight
    
    nr_time_steps = (tend-t0)/dt+1;
    
    nr_v_post_c   = length(V_post_c_vec);
    w_final       = zeros(nr_v_post_c, 1);
    
    for V_post_c_idx = 1:nr_v_post_c
        % Initialization
        V_post_c = V_post_c_vec(V_post_c_idx);
        V_post = V_post_c * ones(nr_time_steps, 1); % Postsyn. membrane pot.
        NMDAr_frac = zeros(nr_time_steps, 1);       % Open NMDAr frac
        NMDAr_cal_cur = zeros(nr_time_steps, 1);    % NMDAr calcium current
        Ca = zeros(nr_time_steps, 1);               % Calcium level
        w  = w_init*ones(nr_time_steps, 1);         % Syanptic weight
        
        % Initialize presynaptic spike times
        t_pre_spikes = calculate_spike_times(...
            pre_spike_freq, nr_pre_spikes, [t0 tend]);
        % Initialize number of "closed NMDArs" right before each spike
        closed_NMDAr_frac_before_spikes = ones(nr_pre_spikes, 1);

        % Model
        tSIM = tic;
        last_pre_spike_time = inf;
        for i = 2:nr_time_steps
            t_next = t0+i*dt;
            
            % Find the most recent presynaptic spike time indices
            recent_t_pre_spikes = t_pre_spikes( t_pre_spikes < t_next );
            spike_number = length(recent_t_pre_spikes);
            % Record number of closed NMDArs if a new spike is coming
            if recent_t_pre_spikes(end) ~= last_pre_spike_time
                closed_NMDAr_frac_before_spikes(spike_number) = ...
                    1-NMDAr_frac(i-1);
            end
            
            % Update variables
            [NMDAr_cal_cur(i), NMDAr_frac(i)] = NMDAr_calcium_current( ...
                t_next, recent_t_pre_spikes, V_post(i), ...
                closed_NMDAr_frac_before_spikes(1:spike_number), ...
                'NMDAr_I_f', NMDAr_I_f);
            Ca(i) = update_Ca(Ca(i-1), NMDAr_cal_cur(i), dt);
            w(i)  = update_w(w(i-1), Ca(i), dt);
            
            % Record last presynaptic spike time
            last_pre_spike_time = recent_t_pre_spikes(end);
        end
        
        % Plot results
        if nr_v_post_c == 1
            tt = linspace(t0, tend, nr_time_steps)';
            if plot_only_w
                plot(tt, w);
                h_leg = legend('Synaptic weight', 'Location', 'Best');
                set(h_leg, 'FontSize', 15);
                xlabel('Time (ms)', 'FontSize', 15);
            else
                figure,
                subplot(3,1,1),
                plot(tt, Ca);
                ylim([0 max([1 max(Ca)])]);
                h_leg = legend('Calcium level', 'Location', 'Best');
                set(h_leg, 'FontSize', 15);
                xlabel('Time (ms)', 'FontSize', 15);
                hold on;
                plot([t0 tend], [0.35 0.35], 'k--');
                plot([t0 tend], [0.55 0.55], 'k--');
                
                subplot(3,1,2),
                plot(tt, w);
                h_leg = legend('Synaptic weight', 'Location', 'Best');
                set(h_leg, 'FontSize', 15);
                xlabel('Time (ms)', 'FontSize', 15);
                
                subplot(3,1,3),
                plot(tt, NMDAr_cal_cur);
            end
        end
        
        % Calculate final weight as the average of min and max in the last
        % 1 percent of simulation
        w_final(V_post_c_idx) = ...
            ( min(w(end-(nr_time_steps-1)/100:end)) + ...
            max(w(end-(nr_time_steps-1)/100:end)) ) / 2;
        fprintf('Post. syn. potential = %i,\t', V_post_c);
        fprintf('w_final/w_init = %1.1f\n', w_final(V_post_c_idx)/w_init);
        toc(tSIM);
    end
    
    % Plot potential vs w(final)/w(0) at the end
    if nr_v_post_c > 1
        figure,
        plot(V_post_c_vec, w_final/w_init);
        hold on;
        plot([-70 -10], [1 1], 'k--');
        xlabel('mV', 'FontSize', 15);
        xlim([-70 -10]);
        ylabel('w(final)/w(init)', 'FontSize', 15);
        ylim([0 4]);
        title('pairing', 'FontSize', 15);
    end
    
    varargout{1} = w_final/w_init;
end

%% Varying the Rate of Presynaptic Stimulation
if strcmp(sim_protocol, 'VtRoPS')
    % Reference: Figure 3B in Shouval 2002
    
    % Input arguments
    pre_spike_freq = varargin{2};       % Presynaptic freq (in Hz)
    nr_pre_spikes  = 100;               % # of presynaptic pulses
    NMDAr_I_f     = 0.50;              % NMDAr fast decay component 
                                        % relative magnitude (0.5)
    w_init         = 0.26;              % Initial weight (0.25-0.26)
    dt             = 0.1;               % in ms
    % --------------------------------------------------------------------
    i = 3;
    while i<=length(varargin),
        switch varargin{i},
            case 'NMDAr_I_f',              NMDAr_I_f = varargin{i+1};
            case 'w_init',                    w_init = varargin{i+1};
            case 'dt',                            dt = varargin{i+1};
            otherwise,
                display(varargin{i});
                error('Unexpected inputs!!!');
        end
        i = i+2;
    end
    % --------------------------------------------------------------------     
    t0             = 0;                 % in ms
    V_rest         = -65;               % Resting potential
    EPSP_norm      = 0.6968;            % See Shouval 2002 supplementary
    EPSP_param_set = 2;
    
    nr_freq = length(pre_spike_freq);
    w_final = zeros(nr_freq,1);
    for freq_idx = 1:nr_freq
        % Calculate spike times, update tend if necessary
        [t_pre_spikes, tend] = calculate_spike_times(...
            pre_spike_freq(freq_idx), nr_pre_spikes, [t0 0]);
        % Initialize number of "closed NMDArs" right before each spike
        closed_NMDAr_frac_before_spikes = ones(nr_pre_spikes, 1);

        % Initialization
        nr_time_steps = (tend-t0)/dt+1;
        V_post = V_rest * ones(nr_time_steps, 1);   % Postsyn. membrane pot.
        NMDAr_frac = zeros(nr_time_steps, 1);       % Open NMDAr frac
        NMDAr_cal_cur = zeros(nr_time_steps, 1);    % NMDAr calcium current
        Ca = zeros(nr_time_steps, 1);               % Calcium level
        w  = w_init*ones(nr_time_steps, 1);         % Syanptic weight
        
        % Model
        tSIM = tic;                 % Keep track of the whole simulations
        tDisplay = tic;             % Keep track display renewal
        display_time_interval = 60; % seconds
        last_pre_spike_time = [];
        for i = 2:nr_time_steps
            t_next = t0+i*dt;
            % Find the most recent presynaptic spike time indices
            recent_t_pre_spikes = t_pre_spikes( t_pre_spikes < t_next );
            spike_number = length(recent_t_pre_spikes);
            % Record number of closed NMDArs if a new spike is coming
            if recent_t_pre_spikes(end) ~= last_pre_spike_time
                closed_NMDAr_frac_before_spikes(spike_number) = ...
                    1-NMDAr_frac(i-1);
            end
            V_post(i) = V_rest + ...
                EPSP(t_next, recent_t_pre_spikes, ...
                EPSP_norm, EPSP_param_set);
            [NMDAr_cal_cur(i), NMDAr_frac(i)] = ...
                NMDAr_calcium_current( ...
                t_next, recent_t_pre_spikes, V_post(i), ...
                closed_NMDAr_frac_before_spikes(1:spike_number), ...
                'NMDAr_I_f', NMDAr_I_f);
            Ca(i) = update_Ca(Ca(i-1), NMDAr_cal_cur(i), dt);
            w(i)  = update_w(w(i-1), Ca(i), dt);
            
            % Display simulation progress
            if toc(tDisplay) > display_time_interval
                % Renew display clock
                tDisplay = tic;
                fprintf('%3.2f%% is completed in %3.2f minutes...\n', ...
                    100*i/nr_time_steps, toc(tSIM)/60);
            end
                        
            % Record last presynaptic spike time
            last_pre_spike_time = recent_t_pre_spikes(end);
        end
        
        % Plot results
        if nr_freq==1
            tt = linspace(t0, tend, nr_time_steps)';
            plot_only_w = 0;
            if plot_only_w
          Figure 5A       plot(tt, w);
                h_leg = legend('Synaptic weight', 'Location', 'Best');
                set(h_leg, 'FontSize', 15);
                xlabel('Time (ms)', 'FontSize', 15);
            else
                figure,
                subplot(2,2,1),
                plot(tt, Ca);
                ylim([0 max([1 max(Ca)])]);
                title('Calcium level', 'FontSize', 15);
                xlabel('Time (ms)', 'FontSize', 15);
                hold on;
                plot([t0 tend], [0.35 0.35], 'k--');
                plot([t0 tend], [0.55 0.55], 'k--');
                
                subplot(2,2,2),
                plot(tt, w);
                title('Synaptic weight', 'FontSize', 15);
                xlabel('Time (ms)', 'FontSize', 15);
                
                subplot(2,2,3),
                plot(tt, NMDAr_cal_cur);
                title('NMDAr calcium current', 'FontSize', 15);
                xlabel('Time (ms)', 'FontSize', 15);
                
                subplot(2,2,4),
                plot(tt, V_post);
                title('Postsynaptic potential', 'FontSize', 15);
                xlabel('Time (ms)', 'FontSize', 15);
                
                figure,
                plot(tt, NMDAr_frac);
                hold on;
                plot(t_pre_spikes, closed_NMDAr_frac_before_spikes, '*');                
                legend('open NMDAr fraction', ...
                    'closed NMDAr fraction right before spikes', ...
                    'FontSize', 15);
                xlabel('Time (ms)', 'FontSize', 15);
            end
        end
        
        % Calculate final weight as the average of min and max in the last
        % 1 percent of simulation
        one_percent = round((nr_time_steps-1)/100);
        w_final(freq_idx) = ( min(w(end-one_percent:end)) + ...
            max(w(end-one_percent:end)) ) / 2;
        fprintf('Pre. stimulation rate = %i,\t', ...
            pre_spike_freq(freq_idx));
        fprintf('w_final/w_init = %1.1f\n', w_final(freq_idx)/w_init);
    end
    
    % Plot presynaptic freq vs w(final)/w(0) at the end
    if nr_freq > 1
        figure,
        plot(pre_spike_freq, w_final/w_init);
        title('frequency', 'FontSize', 15);
    end
    
    % Output arguments
    varargout{1} = w_final/w_init;
end

%% Effect of Varying NMDAR-Binding Kinetics.
if strcmp(sim_protocol, 'EoVNBK-freq')
    % Reference: Figure 5B in Shouval 2002
    % Figure 5B: CaDP('EoVNBK-freq', [0.5:20], 100);
    
    % Input arguments
    pre_spike_freq = varargin{2};       % Presynaptic freq (in Hz)
    nr_pre_spikes  = varargin{3};       % # of presynaptic pulses    
    option = 1;
    
    % Initialization
    w_final_over_w_init_matrix = zeros(length(pre_spike_freq), 1);
    
    for i = 1:3
        if option == 1
            NMDAr_I_f = 0.25*i;
            w_init = 0.26+(2-i)*0.005;
        elseif option == 2
            NMDAr_I_f = 0.50;
            w_init = 0.26+(2-i)*0.005;
        end
        w_final_over_w_init_matrix(:,i) = CaDP('VtRoPS', ...
            pre_spike_freq, nr_pre_spikes, NMDAr_I_f, w_init);
    end
    
    % Plot results
    figure,
    plot(pre_spike_freq, w_final_over_w_init_matrix(:,1), 'r');
    hold on;
    plot(pre_spike_freq, w_final_over_w_init_matrix(:,2), 'k--');
    plot(pre_spike_freq, w_final_over_w_init_matrix(:,3), 'g');   
    plot([0 20], [1 1], 'k--');
    title('frequency', 'FontSize', 15);
    xlabel('Hz', 'FontSize', 15);
    
elseif strcmp(sim_protocol, 'EoVNBK-pairing')
    % Reference: Figure 5A in Shouval 2002
    % Figure 5A: CaDP('EoVNBK-pairing');
    
    % Input arguments
    pre_spike_freq = 1;                 % Presynaptic freq (in Hz)
    nr_pre_spikes  = 100;               % # of presynaptic pulses    
    V_post_c_vec   = (-65:1:-10);       % Fixed post. potential (in mV)
    option         = 1;
    
    % Initialization
    w_final_over_w_init_matrix = zeros(length(V_post_c_vec), 1);
    
    % Model
    for i = 1:3
        if option == 1
            NMDAr_I_f = 0.25*i;
            w_init = 0.25;
        elseif option == 2
            NMDAr_I_f = 0.50;
            w_init = 0.26+(2-i)*0.005;
        end
        w_final_over_w_init_matrix(:,i) = CaDP('PPSwPVC', ...
            [0 1e5], [pre_spike_freq nr_pre_spikes], ...
            V_post_c_vec, NMDAr_I_f, w_init);
    end
    
    % Plot results
    figure,
    plot(V_post_c_vec, w_final_over_w_init_matrix(:,1), 'r');
    hold on;
    plot(V_post_c_vec, w_final_over_w_init_matrix(:,2), 'k--');
    plot(V_post_c_vec, w_final_over_w_init_matrix(:,3), 'g');
    plot([-70 -10], [1 1], 'k--');
    xlabel('mV', 'FontSize', 15);
    xlim([-70 -10]);
    ylabel('w(final)/w(init)', 'FontSize', 15);
    ylim([0 4]);
    title('pairing', 'FontSize', 15);
    
elseif strcmp(sim_protocol, 'EoVNBK-STDP')
    % Reference: Figure 5C in Shouval 2002    
    % Figure 5C: CaDP('EoVNBK-STDP', [-150:5:-51, -50:150, 155:5:250]);
    
    % Input arguments
    spike_timing_diff = varargin{2}; % Spike timing diff. vector (in ms)    
    
    % 0 will be removed from the spike timing diff vector
    spike_timing_diff(spike_timing_diff==0) = [];
    nr_spike_timing_diff = length(spike_timing_diff);
    
    % Initialization
    w_final_over_w_init_matrix = zeros(nr_spike_timing_diff, 3);
    
    % Model
    for i = 1:3
        NMDAr_I_f = 0.25*i;
        %        w_init = 0.25;
        w_final_over_w_init_matrix(:,i) = CaDP('VST', ...
            spike_timing_diff, 'NMDAr_I_f', NMDAr_I_f);
    end
    
    % Plot results
    figure,
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,1), 'r');
    hold on;
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,2), 'k--');
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,3), 'g');
    plot([min(spike_timing_diff) max(spike_timing_diff)], [1 1], 'k:');
    xlabel('\Deltat (ms)', 'FontSize', 15);
    xlim([min(spike_timing_diff) max(spike_timing_diff)]);
    ylabel('w(final)/w(init)', 'FontSize', 15);
    ylim([0.2 2.3]);
    title('STDP', 'FontSize', 15);

end

%% Induction of Synaptic Plasticity by Varying Spike Timing
if strcmp(sim_protocol, 'VST')
    % Reference: Figure 3C in Shouval 2002
    % Figure 3C: CaDP('VST', [-150:200], 0.50);   
    % Figure S9: CaDP('VST', [-120:10:-20 -20:3:29 30:10:100], 'dt', 0.5,
    %                 'NMDAr_I_f', 1, 'BPAP_tau_s', 75);

    % Input arguments
    spike_timing_difference = varargin{2}; % t_post - t_pre
    NMDAr_I_f               = 0.50;        % NMDAr fast decay magnitude
    NMDAr_tau_s             = 200;         % NMDAr slow decay rate
    BPAP_tau_s              = 25;          % BPAP slow decay rate
    pre_spike_freq          = 1;           % Presynaptic freq (in Hz)
    nr_pre_spikes           = 100;         % # of presynaptic pulses
    w_init                  = 0.25;        % Initial weight
    dt                      = 0.1;         % in ms
    % --------------------------------------------------------------------
    i = 3;
    while i<=length(varargin),
        switch varargin{i},
            case 'NMDAr_I_f',              NMDAr_I_f = varargin{i+1};
            case 'NMDAr_tau_s',          NMDAr_tau_s = varargin{i+1};
            case 'BPAP_tau_s',            BPAP_tau_s = varargin{i+1};
            case 'pre_spike_freq',    pre_spike_freq = varargin{i+1};
            case 'w_init',                    w_init = varargin{i+1};
            case 'dt',                            dt = varargin{i+1};
            case 'nr_pre_spikes',      nr_pre_spikes = varargin{i+1};
            otherwise,
                display(varargin{i});
                error('Unexpected inputs!!!');
        end
        i = i+2;
    end
    % -------------------------------------------------------------------- 
    t0             = 0;                 % in ms
    tend           = 1;                 % in ms (will be changed)
    V_rest         = -65;               % Resting potential
    EPSP_norm      = 0.6968;            % See Shouval 2002 supplementary
    EPSP_param_set = 1;
    BPAP_type      = 'BPAP + ADP';
    plot_only_w    = 0;
                
    % Initialize presynaptic spike times
    [t_pre_spikes, tend] = calculate_spike_times(...
        pre_spike_freq, nr_pre_spikes, [t0 tend]);
    % Initialize number of "closed NMDArs" right before each spike
    closed_NMDAr_frac_before_spikes = ones(nr_pre_spikes, 1);
    
    % Remove 0 from spike timing difference vector
    if ismember(0, spike_timing_difference)
        spike_timing_difference = setdiff(spike_timing_difference, 0);
    end
    
    % Find the length of the spike timing difference vector
    nr_timing_diffs = length(spike_timing_difference);
    % Initialize final weight vector
    w_final = zeros(nr_timing_diffs, 1);
    for timing_diff_idx = 1:nr_timing_diffs
        % Find the current spike timing difference
        spt = spike_timing_difference(timing_diff_idx);
        
        % Initialize postsynaptic spike times
        if spt > 0
            t_post_spikes = t_pre_spikes + abs(spt);
        elseif spt < 0
            t_post_spikes = t_pre_spikes;
            t_pre_spikes  = t_pre_spikes + abs(spt);
        end
        tend = tend + spt;
        
        nr_time_steps = (tend-t0)/dt+1;
        
        % Initialize
        V_post = V_rest * ones(nr_time_steps, 1);   % Postsyn. membrane pot.
        NMDAr_frac = zeros(nr_time_steps, 1);       % Open NMDAr frac
        NMDAr_cal_cur = zeros(nr_time_steps, 1);    % NMDAr calcium current
        Ca = zeros(nr_time_steps, 1);               % Calcium level
        w  = w_init*ones(nr_time_steps, 1);         % Syanptic weight
        
        % Model
        tSIM = tic;                 % Keep track of the whole simulations
        tDisplay = tic;             % Keep track display renewal
        display_time_interval = 6;  % seconds
        last_pre_spike_time = inf;
        for i = 2:nr_time_steps
            t_next = t0+i*dt;
            
            % Find the most recent presynaptic spike time indices
            recent_t_pre_spikes = t_pre_spikes( t_pre_spikes < t_next );
            spike_number = length(recent_t_pre_spikes);
            % Record number of closed NMDArs if a new spike is coming
            if ~isempty(recent_t_pre_spikes) && ...
                    recent_t_pre_spikes(end) ~= last_pre_spike_time
                closed_NMDAr_frac_before_spikes(spike_number) = ...
                    1-NMDAr_frac(i-1);
            end
            % Record last presynaptic spike time
            if isempty(recent_t_pre_spikes)
                last_pre_spike_time = inf;
            else
                last_pre_spike_time = recent_t_pre_spikes(end);
            end
            
            % Calculate potential at the dentrite
            V_post(i) = V_rest + ...
                EPSP(t_next, t_pre_spikes, EPSP_norm, ...
                EPSP_param_set) + ...
                BPAP(t_next, t_post_spikes, 'BPAP_type', BPAP_type, ...
                'BPAP_tau_s', BPAP_tau_s);
            % Calculate NMDA calcium current
            [NMDAr_cal_cur(i), NMDAr_frac(i)] = ...
                NMDAr_calcium_current( ...
                t_next, t_pre_spikes, V_post(i), ...
                closed_NMDAr_frac_before_spikes(1:spike_number), ...
                'NMDAr_tau_s', NMDAr_tau_s, 'NMDAr_I_f', NMDAr_I_f);
            % Update Calcium level
            Ca(i) = update_Ca(Ca(i-1), NMDAr_cal_cur(i), dt);
            % Update syanptic weight
            w(i)  = update_w(w(i-1), Ca(i), dt);
            
            % Display simulation progress
            
            if toc(tDisplay) > display_time_interval
                % renew display clock
                tDisplay = tic;
                % display progress
                fprintf('%2.2f%% is completed in %3.2f minutes...\n', ...
                    100*i/nr_time_steps, toc(tSIM)/60);
            end
        end
        
        % Calculate final weight as the average of min and max in the last
        % 1 percent of simulation
        one_percent = round((nr_time_steps-1)/100);
        w_final(timing_diff_idx) = ( min(w(end-one_percent:end)) + ...
            max(w(end-one_percent:end)) ) / 2;
        fprintf('Spike timing difference = %i,\t', spt);
        fprintf('w_final/w_init = %1.1f\n', ...
            w_final(timing_diff_idx)/w_init);
        toc(tSIM);

    end
        
    % Plot results
    if nr_timing_diffs == 1
        tt = linspace(t0, tend, nr_time_steps)';
        if plot_only_w
            plot(tt, w);
            h_leg = legend('Synaptic weight', 'Location', 'Best');
            set(h_leg, 'FontSize', 15);
            xlabel('Time (ms)', 'FontSize', 15);
        else
            figure,
            subplot(4,1,1),
            plot(tt, Ca);
            ylim([0 max([1 max(Ca)])]);
            h_leg = legend('Calcium level', 'Location', 'Best');
            set(h_leg, 'FontSize', 15);
            xlabel('Time (ms)', 'FontSize', 15);
            hold on;
            plot([t0 tend], [0.35 0.35], 'k--');
            plot([t0 tend], [0.55 0.55], 'k--');
            
            subplot(4,1,2),
            plot(tt, w);
            h_leg = legend('Synaptic weight', 'Location', 'Best');
            set(h_leg, 'FontSize', 15);
            
            subplot(4,1,3),
            plot(tt, NMDAr_cal_cur);
            h_leg = legend('NMDAr calcium current', 'Location', 'Best');
            set(h_leg, 'FontSize', 15);
            
            subplot(4,1,4),
            plot(tt, V_post);
            h_leg = legend('Post potential', 'Location', 'Best');
            set(h_leg, 'FontSize', 15);            
            
            figure,
            AX = plotyy(tt, NMDAr_frac, tt, V_post);
            h_leg = legend('NMDAr fraction', 'Post potential', ...
                'Location', 'Best');
            set(h_leg, 'FontSize', 15);
        end
    end
    
    % Plot potential vs w(final)/w(0) at the end
    if nr_timing_diffs > 1
        figure,
        plot(spike_timing_difference, w_final/w_init);
        xlabel('\Delta t = t_p_o_s_t - t_p_r_e', 'FontSize', 15);
        ylabel('w(final)/w(init)', 'FontSize', 15);
        title('STDP', 'FontSize', 15);
        %         hold on;
        %         plot([-70 -10], [1 1], 'k--');
        %         xlim([-70 -10]);
        %         ylim([0 4]);
    end
    
    % Output arguments
    varargout{1} = w_final/w_init;
        
end

%% Effect of BPAP on STDP
if strcmp(sim_protocol, 'EoBPAPoSTDP')
    % Reference: Figure 4A in Shouval 2002
    % Figure 4A: EoBPAPoSTDP_dt_0pnt5_Mar30_0224.mat
    
    % Input arguments
    spike_timing_diff = varargin{2}; % Spike timing diff. vector (in ms)    
    BPAP_tau_s_vec    = [15 25 50];  % BPAP slow component decay rate
    dt                = 0.1;
    save_workspace    = 0;
    % --------------------------------------------------------------------
    i = 3;
    while i<=length(varargin),
        switch varargin{i},
            case 'dt',                            dt = varargin{i+1};
            case 'save_workspace',    save_workspace = varargin{i+1};
            otherwise,
                display(varargin{i});
                error('Unexpected inputs!!!');
        end
        i = i+2;
    end
    % --------------------------------------------------------------------
    % 0 will be removed from the spike timing diff vector
    spike_timing_diff(spike_timing_diff==0) = [];
    nr_spike_timing_diff = length(spike_timing_diff);
    
    % Initialization
    w_final_over_w_init_matrix = zeros(nr_spike_timing_diff, 3);
    
    % Model
    for i = 1:3
        BPAP_tau_s = BPAP_tau_s_vec(i);
        w_final_over_w_init_matrix(:,i) = CaDP('VST', ...
            spike_timing_diff, 'BPAP_tau_s', BPAP_tau_s, 'dt', dt);
    end
    
    % Save_results
    if save_workspace
        foldername = 'saved_workspaces/';
        filename = [sim_protocol datestr(now,'_mmmdd_HHMM') '.mat'];
        save([foldername filename]);
    end
    
    % Plot results
    figure,
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,1), 'r');
    hold on;
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,2), 'k--');
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,3), 'g');
    plot([min(spike_timing_diff) max(spike_timing_diff)], [1 1], 'k:');
    xlabel('\Deltat (ms)', 'FontSize', 15);
    xlim([min(spike_timing_diff) max(spike_timing_diff)]);
    ylabel('w(final)/w(init)', 'FontSize', 15);
    ylim([0.2 2.3]);
    title('STDP', 'FontSize', 15);
    
end

%% Effect of stimulation rate on STDP
if strcmp(sim_protocol, 'EoSRoSTDP')
    % Reference: Figure 4B in Shouval 2002
    % Figure 4B:
    
    % Input arguments
    spike_timing_diff  = varargin{2};   % Spike timing diff. vector (in ms)
    pre_spike_freq_vec = [1 5 10];      % Presynaptic freq vec (in Hz)
    w_init             = [.25 .25 .26]; % Initial weight
    dt                 = 1;
    save_workspace     = 0;
    % --------------------------------------------------------------------
    i = 3;
    while i<=length(varargin),
        switch varargin{i},
            case 'dt',                            dt = varargin{i+1};
            case 'save_workspace',    save_workspace = varargin{i+1};
            otherwise,
                display(varargin{i});
                error('Unexpected inputs!!!');
        end
        i = i+2;
    end
    % --------------------------------------------------------------------

    % 0 will be removed from the spike timing diff vector
    spike_timing_diff(spike_timing_diff==0) = [];
    nr_spike_timing_diff = length(spike_timing_diff);
    
    % Initialization
    w_final_over_w_init_matrix = zeros(nr_spike_timing_diff, 3);
    
    % Model
    for i = 1:3
        w_final_over_w_init_matrix(:,i) = CaDP('VST', ...
            spike_timing_diff, ...
            'pre_spike_freq', pre_spike_freq_vec(i), ...
            'w_init', w_init(i), 'dt', dt);
    end
    
    % Save_results
    if save_workspace
        foldername = 'saved_workspaces/';
        filename = [sim_protocol datestr(now,'_mmmdd_HHMM') '.mat'];
        save([foldername filename]);
    end

    % Plot results
    figure,
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,1), 'b');
    hold on;
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,2), 'r');
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,3), 'k');
    plot([min(spike_timing_diff) max(spike_timing_diff)], [1 1], 'k:');
    plot([0 0], [0 4], 'k:');
    xlabel('\Deltat (ms)', 'FontSize', 15);
    xlim([min(spike_timing_diff) max(spike_timing_diff)]);
    ylabel('w(final)/w(init)', 'FontSize', 15);
    ylim([0 4.2]);
    title('STDP', 'FontSize', 15);
    
end

%% Effect of NMDAR slow decay rate
if strcmp(sim_protocol, 'EoNMDARSDR')
    % Reference: Not from Shouval 2002
    % EoNMDARSDR_dt_0pnt5_Mar30_0712.mat
    
    % Input arguments
    spike_timing_diff  = varargin{2};   % Spike timing diff. vector (in ms)
    NMDAr_tau_s_vec    = [100 200 300]; % NMDAR slow decay rate vector
    dt                 = 0.1;
    save_workspace     = 0;
    % --------------------------------------------------------------------
    i = 3;
    while i<=length(varargin),
        switch varargin{i},
            case 'dt',                            dt = varargin{i+1};
            case 'save_workspace',    save_workspace = varargin{i+1};
            otherwise,
                display(varargin{i});
                error('Unexpected inputs!!!');
        end
        i = i+2;
    end
    % --------------------------------------------------------------------

    % 0 will be removed from the spike timing diff vector
    spike_timing_diff(spike_timing_diff==0) = [];
    nr_spike_timing_diff = length(spike_timing_diff);
    
    % Initialization
    w_final_over_w_init_matrix = zeros(nr_spike_timing_diff, 3);
    
    % Model
    for i = 1:3
        w_final_over_w_init_matrix(:,i) = CaDP('VST', ...
            spike_timing_diff, 'NMDAr_tau_s', NMDAr_tau_s_vec(i), ...
            'dt', dt);
    end
    
    % Save_results
    if save_workspace
        foldername = 'saved_workspaces/';
        filename = [sim_protocol datestr(now,'_mmmdd_HHMM') '.mat'];
        save([foldername filename]);
    end

    % Plot results
    figure,
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,1), 'r');
    hold on;
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,2), 'k--');
    plot(spike_timing_diff, w_final_over_w_init_matrix(:,3), 'g');
    for i = 3:-1:1
        plot_leg{i} = ['NMDAr \tau_s=' num2str(NMDAr_tau_s_vec(i))];
    end
    h_leg = legend(plot_leg{1}, plot_leg{2}, plot_leg{3}, ...
        'Location', 'Best');
    set(h_leg, 'FontSize', 15);
    plot([min(spike_timing_diff) max(spike_timing_diff)], [1 1], 'k:');
    plot([0 0], [0.2 2.2], 'k:');
    xlabel('\Deltat (ms)', 'FontSize', 15);
    xlim([min(spike_timing_diff) max(spike_timing_diff)]);
    ylabel('w(final)/w(init)', 'FontSize', 15);
    ylim([0.2 2.2]);
    title('STDP', 'FontSize', 15);
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% STOCHASTIC CADP, Shouval 2004 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mean Calcium transient and their standard deviation
if strcmp(sim_protocol, 'sCDAP_MuSigmaCa')
    % Reference: Figure 1C in Shouval 2004
    
    % Input arguments
    spike_timing_diff = varargin{2};    % t_post - t_pre
    NMDAr_version     = varargin{3};    % 'deterministic-with-stochastic-parameters'
                                        % 'stochastic'
    nr_sims           = varargin{4};    % # of simulations
    sigma_square      = varargin{5};    % variance of G_NMDAr
    
    dt             = 0.1;         % in ms
    BPAP_delay     = 0;
    Z              = 10;          % # on NMDA receptors
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
        G_NMDA_for_each_spike          = ...
            generate_random_maximal_chord_conductance(...
            spike_timing_diff, Z, 'sigma_square', sigma_square);
        
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
                'G_NMDA_for_each_spike', G_NMDA_for_each_spike);
            % Update Calcium level
            Ca(i) = update_Ca(Ca(i-1), NMDAr_cal_cur(i), dt);
        end
        
        % Record Calcium peak
        Ca_peaks(k) = max(Ca);
        
        % Display simulation progress
        tDISP = display_sim_progress(k, nr_sims, tDISP, tSIM);
    end
    
    % Plot results
    figure,
    if nr_sims == 1        
        tt = linspace(t0, tend, nr_time_steps)';
        plot(tt, Ca);
        h_leg = legend('Calcium level', 'Location', 'Best');
        set(h_leg, 'FontSize', 15);
        xlabel('Time (ms)', 'FontSize', 15);
    elseif nr_sims > 1
        histogram(Ca_peaks, 'Normalization', 'probability');
        CV_sim = std(nonzeros(Ca_peaks))/mean(nonzeros(Ca_peaks));
        CV_theory = 0.095 + ...
            0.0045*max([0 spike_timing_diff]) - ...
            0.00067*min([0 spike_timing_diff]);
        CV_error = CV_sim - CV_theory;
        title(['\Deltat = ' num2str(spike_timing_diff) ...
            ', nonzero Ca peaks CV = ', ...
            num2str(CV_sim, '%1.2f'), ...
            ', sim-theory =  ', num2str(CV_error, '%1.2f')]);
    end
    
end

%% Deterministic/Stochastic STDP with Shouval 2004 parameters
if strcmp(sim_protocol, 'sCDAP_VST')
    % Reference: Figure 1B in Shouval 2004

    % Input arguments
    spike_timing_difference = varargin{2}; % t_post - t_pre
    NMDAr_version           = varargin{3}; 
    % 'deterministic-with-stochastic-parameters'
    % 'stochastic'
    plot_results   = 1;
    pre_spike_freq = 1;           % Presynaptic freq (in Hz)
    nr_pre_spikes  = 100;         % # of presynaptic pulses
    w_init         = 0.25;        % Initial weight
    dt             = 0.1;         % in ms
    BPAP_delay     = 0;
    Z              = 10;          % # on NMDA receptors
    % --------------------------------------------------------------------
    i = 4;
    while i<=length(varargin),
        switch varargin{i},
            case 'sigma_square',        sigma_square = varargin{i+1};
            case 'plot_results',        plot_results = varargin{i+1};
            case 'Z',                              Z = varargin{i+1};
            case 'pre_spike_freq',    pre_spike_freq = varargin{i+1};
            case 'w_init',                    w_init = varargin{i+1};
            case 'dt',                            dt = varargin{i+1};
            case 'nr_pre_spikes',      nr_pre_spikes = varargin{i+1};
            otherwise,
                display(varargin{i});
                error('Unexpected inputs!!!');
        end
        i = i+2;
    end
    % -------------------------------------------------------------------- 
    t0             = 0;                 % in ms
    tend           = 1;                 % in ms (will be changed)
    V_rest         = -65;               % Resting potential
    EPSP_norm      = 0.6968;            % See Shouval 2002 supplementary
    EPSP_param_set = 1;
    plot_only_w    = 0;
                
    % Initialize presynaptic spike times
    [t_pre_spikes, tend] = calculate_spike_times(...
        pre_spike_freq, nr_pre_spikes, [t0 tend]);
    % Initialize number of "closed NMDArs" right before each spike
    closed_NMDAr_frac_before_spikes = ones(nr_pre_spikes, 1);
    G_NMDA_for_each_spike = zeros(nr_pre_spikes, 1);
    
    % Remove 0 from spike timing difference vector
    if ismember(0, spike_timing_difference)
        spike_timing_difference = setdiff(spike_timing_difference, 0);
    end
    
    % Find the length of the spike timing difference vector
    nr_timing_diffs = length(spike_timing_difference);
    % Initialize final weight vector
    w_final = zeros(nr_timing_diffs, 1);
    for timing_diff_idx = 1:nr_timing_diffs
        % Find the current spike timing difference
        spt = spike_timing_difference(timing_diff_idx);
        
        % Initialize postsynaptic spike times
        if spt > 0
            t_post_spikes = t_pre_spikes + abs(spt);
        elseif spt < 0
            t_post_spikes = t_pre_spikes;
            t_pre_spikes  = t_pre_spikes + abs(spt);
        end
        tend = tend + spt;
        
        nr_time_steps = (tend-t0)/dt+1;
        
        % Initialize
        V_post = V_rest * ones(nr_time_steps, 1);   % Postsyn. membrane pot.
        NMDAr_frac = zeros(nr_time_steps, 1);       % Open NMDAr frac
        NMDAr_cal_cur = zeros(nr_time_steps, 1);    % NMDAr calcium current
        Ca = zeros(nr_time_steps, 1);               % Calcium level
        w  = w_init*ones(nr_time_steps, 1);         % Syanptic weight
        
        % Model
        tSIM = tic;                 % Keep track of the whole simulations
        tDisplay = tic;             % Keep track display renewal
        display_time_interval = 6;  % seconds
        last_pre_spike_time = inf;
        for i = 2:nr_time_steps
            t_next = t0+i*dt;
            
            % Find the most recent presynaptic spike time indices
            recent_t_pre_spikes = t_pre_spikes( t_pre_spikes < t_next );
            spike_number = length(recent_t_pre_spikes);
            % Record number of closed NMDArs if a new spike is coming
            if ~isempty(recent_t_pre_spikes) && ...
                    recent_t_pre_spikes(end) ~= last_pre_spike_time
                closed_NMDAr_frac_before_spikes(spike_number) = ...
                    1-NMDAr_frac(i-1);
                if strcmp(NMDAr_version, 'stochastic')
                    G_NMDA_for_each_spike(spike_number) = ...
                        generate_random_maximal_chord_conductance(...
                        spt, Z, 'sigma_square', sigma_square);
                end
            end
            % Record last presynaptic spike time
            if isempty(recent_t_pre_spikes)
                last_pre_spike_time = inf;
            else
                last_pre_spike_time = recent_t_pre_spikes(end);
            end
            
            % Calculate potential at the dentrite
            V_post(i) = V_rest + ...
                EPSP(t_next, t_pre_spikes, EPSP_norm, ...
                EPSP_param_set) + ...
                BPAP(t_next, t_post_spikes, 'version', 'stochastic', ...
                'BPAP_delay', BPAP_delay);
            % Calculate NMDA calcium current
            [NMDAr_cal_cur(i), NMDAr_frac(i)] = ...
                NMDAr_calcium_current( ...
                t_next, t_pre_spikes, V_post(i), ...
                closed_NMDAr_frac_before_spikes(1:spike_number), ...
                'version', NMDAr_version, ...
                'G_NMDA_for_each_spike', ...
                G_NMDA_for_each_spike(1:spike_number));
            % Update Calcium level
            Ca(i) = update_Ca(Ca(i-1), NMDAr_cal_cur(i), dt, ...
                'version', 'stochastic');
            % Update syanptic weight
            w(i)  = update_w(w(i-1), Ca(i), dt, 'version', 'stochastic');
            
            % Display simulation progress            
            if toc(tDisplay) > display_time_interval
                % renew display clock
                tDisplay = tic;
                % display progress
                fprintf('%2.2f%% is completed in %3.2f minutes...\n', ...
                    100*i/nr_time_steps, toc(tSIM)/60);
            end
        end
        
        % Calculate final weight as the average of min and max in the last
        % 1 percent of simulation
        one_percent = round((nr_time_steps-1)/100);
        w_final(timing_diff_idx) = ( min(w(end-one_percent:end)) + ...
            max(w(end-one_percent:end)) ) / 2;
        fprintf('Spike timing difference = %i,\t', spt);
        fprintf('w_final/w_init = %1.3f\n', ...
            w_final(timing_diff_idx)/w_init);
        toc(tSIM);

    end
        
    % Plot results
    if nr_timing_diffs == 1 && plot_results
        tt = linspace(t0, tend, nr_time_steps)';
        if plot_only_w
            plot(tt, w);
            h_leg = legend('Synaptic weight', 'Location', 'Best');
            set(h_leg, 'FontSize', 15);
            xlabel('Time (ms)', 'FontSize', 15);
        else
            figure,
            subplot(4,1,1),
            plot(tt, Ca);
            ylim([0 max([1 max(Ca)])]);
            h_leg = legend('Calcium level', 'Location', 'Best');
            set(h_leg, 'FontSize', 15);
            xlabel('Time (ms)', 'FontSize', 15);
            hold on;
            plot([t0 tend], [0.35 0.35], 'k--');
            plot([t0 tend], [0.55 0.55], 'k--');
            
            subplot(4,1,2),
            plot(tt, w);
            h_leg = legend('Synaptic weight', 'Location', 'Best');
            set(h_leg, 'FontSize', 15);
            
            subplot(4,1,3),
            plot(tt, NMDAr_cal_cur);
            h_leg = legend('NMDAr calcium current', 'Location', 'Best');
            set(h_leg, 'FontSize', 15);
            
            subplot(4,1,4),
            plot(tt, V_post);
            h_leg = legend('Post potential', 'Location', 'Best');
            set(h_leg, 'FontSize', 15);            
            
            figure,
            AX = plotyy(tt, NMDAr_frac, tt, V_post);
            h_leg = legend('NMDAr fraction', 'Post potential', ...
                'Location', 'Best');
            set(h_leg, 'FontSize', 15);
        end
    end
    
    % Plot potential vs w(final)/w(0) at the end
    if nr_timing_diffs > 1 && plot_results
        figure,
        plot(spike_timing_difference, w_final/w_init);
        xlabel('\Delta t = t_p_o_s_t - t_p_r_e', 'FontSize', 15);
        ylabel('w(final)/w(init)', 'FontSize', 15);
        title('STDP', 'FontSize', 15);
        hold on;
        plot([min(spike_timing_difference) ...
            max(spike_timing_difference)], [1 1], 'k:');
        
        ylim([0.4 2.5]);
    end
    
    % Output arguments
    varargout{1} = w_final/w_init;
        
end

%% STDP window in Shouval 2004
if strcmp(sim_protocol, 'sCDAP_VST_det_and_sto')
    % Reference: Figure 2 in Shouval 2004
    
    % Input arguments
    spike_timing_difference = varargin{2}; % t_post - t_pre
    nr_reps                 = varargin{3};
    dt                      = 0.1;   
    Z                       = 10;
    save_workspace          = 0;
    plot_individual_results = 0;
    % --------------------------------------------------------------------
    i = 4;
    while i<=length(varargin),
        switch varargin{i},
            case 'save_workspace',    save_workspace = varargin{i+1};
            case 'Z',                              Z = varargin{i+1};
            case 'dt',                            dt = varargin{i+1};
            otherwise,
                display(varargin{i});
                error('Unexpected inputs!!!');
        end
        i = i+2;
    end
    % -------------------------------------------------------------------- 
    % Remove 0 from spike timing difference vector
    if ismember(0, spike_timing_difference)
        spike_timing_difference = setdiff(spike_timing_difference, 0);
    end
    
    % Find the length of the spike timing difference vector
    nr_timing_diffs = length(spike_timing_difference);

    % Initialized
    sigma_square_vec = zeros(nr_timing_diffs, 1);
    deterministic_w_final_over_w_init = zeros(nr_timing_diffs, 1); 
    stochastic_w_final_over_w_init = zeros(nr_timing_diffs, nr_reps); 
    
    for i = 1:nr_timing_diffs
        % Current spike timing difference 
        spt = spike_timing_difference(i);
        
        % Find variance of G_NMDA gamma distribution variance that matches
        % the theory
        sigma_square_vec(i) = find_sCDAP_G_NMDA_gamma_dist_variance(spt);
        
        % Calculate deterministic synaptic plasticity
        deterministic_w_final_over_w_init(i) = CaDP('sCDAP_VST', spt, ...
            'deterministic-with-stochastic-parameters', ...
            'plot_results', plot_individual_results, ...
            'sigma_square', sigma_square_vec(i), 'dt', dt, 'Z', Z);
        
        % Calculate stochastic synaptic plasticity
        for j = 1:nr_reps
            stochastic_w_final_over_w_init(i, j) = ...
                CaDP('sCDAP_VST', spt, 'stochastic', ...
                'plot_results', plot_individual_results, ...
                'sigma_square', sigma_square_vec(i), 'dt', dt, 'Z', Z);
        end
    end
    
    % Save_results
    if save_workspace
        foldername = 'saved_workspaces/';
        filename = [sim_protocol datestr(now,'_mmmdd_HHMM') '.mat'];
        save([foldername filename]);
    end
    
    % Plot results
    figure,
    plot(spike_timing_difference, deterministic_w_final_over_w_init, 'k--');
    hold on;
    plot(spike_timing_difference, stochastic_w_final_over_w_init, 'b+');
    errorbar(spike_timing_difference, ...
        mean(stochastic_w_final_over_w_init, 2), ...
        std(stochastic_w_final_over_w_init, 0, 2), 'g', ...
        'LineWidth', 2);
    plot([min(spike_timing_difference) ...
        max(spike_timing_difference)], [1 1], 'k:');
    xlim([min(spike_timing_difference)-1 max(spike_timing_difference)+1]);
    xlabel('\Delta t = t_p_o_s_t - t_p_r_e', 'FontSize', 15);
    ylabel('w(final)/w(init)', 'FontSize', 15);
    title(['STDP (Z=' num2str(Z) ')'], 'FontSize', 15);
    xlim([-150 150]);
    
    % Exponential fits
    shift = 1;
    idx = find(spike_timing_difference>0, 1);
    x0  = spike_timing_difference(idx);
    y0  = mean(stochastic_w_final_over_w_init(idx,:), 2) - shift;
    pos_exp_fit_tau = 14;
    pos_exp_fit_C   = y0/exp(-x0/pos_exp_fit_tau);
    Xpos = (x0:max(spike_timing_difference))';
    Ypos = pos_exp_fit_C * exp(-Xpos/14) + shift;
    plot(Xpos, Ypos, 'r');
        
end


end