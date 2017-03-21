function CaDP(varargin)
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
    % Reference: Figure 3 in Shouval 2002
    % Figure 3A: CaDP('PPSwPVC', [0 1e5], [1 100], [-65:1:-10]);
    
    % Input arguments
    t0             = varargin{2}(1);    % in ms
    tend           = varargin{2}(2);    % in ms
    pre_spike_freq = varargin{3}(1);    % Presynaptic freq (in Hz)
    nr_pre_spikes  = varargin{3}(2);    % # of presynaptic pulses
    V_post_c_vec   = varargin{4};       % Fixed post. potential (in mV)
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
        NMDAr_cal_cur = zeros(nr_time_steps, 1);    % NMDAr calcium current
        Ca = zeros(nr_time_steps, 1);               % Calcium level
        w  = w_init*ones(nr_time_steps, 1);         % Syanptic weight
        
        t_pre_spikes = calculate_spike_times(...
            pre_spike_freq, nr_pre_spikes, [t0 tend]);        
        
        % Model
        tSIM = tic;
        for i = 2:nr_time_steps
            t_next = t0+i*dt;
            % Find the most recent presynaptic spike time index
            recent_t_pre_spike_idx = find(t_pre_spikes>t_next,1)-1;
            % Choose inf as the pre spike time if no spikes yet, as this will
            % return a 0 calcium current because of the use of the Heaviside
            % function
            if recent_t_pre_spike_idx==0
                recent_t_pre_spike = inf;
            elseif isempty(recent_t_pre_spike_idx)
                recent_t_pre_spike = t_pre_spikes(end);
            else
                recent_t_pre_spike = t_pre_spikes(recent_t_pre_spike_idx);
            end
            NMDAr_cal_cur(i) = NMDAr_calcium_current( ...
                t_next, recent_t_pre_spike, V_post(i));
            Ca(i) = update_Ca(Ca(i-1), NMDAr_cal_cur(i), dt);
            w(i)  = update_w(w(i-1), Ca(i), dt);
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
end


end