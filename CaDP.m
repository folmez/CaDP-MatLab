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
    dt           = 0.1;              % in ms
    V_rest       = -65;              % in mV
    BPAP_delay   = 2;                % in ms
    plot_NMDAr_calcium_current = 0;
    
    nr_time_steps = (tend-t0)/dt+1;
    
    % Initialization
    V_post = V_rest * ones(nr_time_steps, 1);   % Postsyn. membrane pot.
    NMDAr_cal_cur = zeros(nr_time_steps, 1);    % NMDAr calcium current
    NMDAr_frac = zeros(nr_time_steps, 1);       % Open NMDAr fraction
    Ca = zeros(nr_time_steps, 1);               % Calcium level
    
    % Model
    for i = 2:nr_time_steps
        t_next = t0+i*dt;
        % Calculate potential at the dentrite
        V_post(i) = V_rest + EPSP(t_next, t_pre_spike, 1, 1) + ...
            BPAP(t_next, t_post_spike, BPAP_delay, BPAP_type);
        % Calculate NMDA calcium current
        [NMDAr_cal_cur(i), NMDAr_frac(i)] = NMDAr_calcium_current( ...
            t_next, t_pre_spike, V_post(i));
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
    NMDA_r_I_f     = varargin{5};       % NMDAr fast decay component 
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
                NMDA_r_I_f, ...
                closed_NMDAr_frac_before_spikes(1:spike_number));
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
    nr_pre_spikes  = varargin{3};       % # of presynaptic pulses
    NMDA_r_I_f     = varargin{4};       % NMDAr fast decay component 
                                        % relative magnitude (0.5)
    w_init         = varargin{5};       % Initial weight (0.25-0.26)

    t0             = 0;                 % in ms
    dt             = 1;                 % in ms
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
                t_next, recent_t_pre_spikes, V_post(i), NMDA_r_I_f, ...
                closed_NMDAr_frac_before_spikes(1:spike_number));
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
end

%% Induction of Synaptic Plasticity by Varying Spike Timing
if strcmp(sim_protocol, 'VST')
    % Reference: Figure 3C in Shouval 2002
    % Figure 3C: CaDP('VST', [-150:200], 1);
    
    % Input arguments
    spike_timing_difference = varargin{2};      % t_post - t_pre
    pre_spike_freq          = varargin{3};      % Presynaptic freq (in Hz)
    
    t0             = 0;                 % in ms
    tend           = 1;                 % in ms (will be changed)
    nr_pre_spikes  = 100;               % # of presynaptic pulses
    V_rest         = -65;               % Resting potential
    NMDA_r_I_f     = 0.5;               % NMDAr fast decay component 
    w_init         = 0.25;              % Initial weight
    dt             = 1;                 % in ms
    EPSP_norm      = 0.6968;            % See Shouval 2002 supplementary
    EPSP_param_set = 1;
    BPAP_delay     = 2;                 % in ms
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
        tSIM = tic;
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
                EPSP(t_next, last_pre_spike_time, EPSP_norm, EPSP_param_set) + ...
                BPAP(t_next, t_post_spikes, BPAP_delay, BPAP_type);
            % Calculate NMDA calcium current
            [NMDAr_cal_cur(i), NMDAr_frac(i)] = ...
                NMDAr_calcium_current( ...
                t_next, recent_t_pre_spikes, V_post(i), NMDA_r_I_f, ...
                closed_NMDAr_frac_before_spikes(1:spike_number));
            % Update Calcium level
            Ca(i) = update_Ca(Ca(i-1), NMDAr_cal_cur(i), dt);
            % Update syanptic weight
            w(i)  = update_w(w(i-1), Ca(i), dt);
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
            
            subplot(3,1,3),
            plot(tt, NMDAr_cal_cur);
            h_leg = legend('NMDAr calcium current', 'Location', 'Best');
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


end