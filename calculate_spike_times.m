function varargout = calculate_spike_times(varargin)

% Input arguments
pre_spike_freq = varargin{1};    % Presynaptic freq (in Hz)
nr_pre_spikes  = varargin{2};    % # of presynaptic pulses
t0             = varargin{3}(1); % in ms
tend           = varargin{3}(2); % in ms

disable_warnings = 0;
t_unit           = 1000;         % s to ms conversion factor

% Check inputs
if nr_pre_spikes/((tend-t0)/t_unit) > pre_spike_freq
    prev_tend = tend;
    tend = t_unit*nr_pre_spikes/pre_spike_freq + t0;
    tend = ceil(tend);
    if prev_tend ~= 0 && ~disable_warnings
        warning_message = ['FO: Simulation end time update: ' ...
            'Before = ' num2str(prev_tend, '%1.1e') ...
            ', Now = ' num2str(tend,'%1.1e')];
        warning(warning_message);
    end
end
    
% Model
inter_spike_int = t_unit/pre_spike_freq;                       % in ms
t_pre_spike_times = t0 + (0:nr_pre_spikes-1)'*inter_spike_int; % in ms

% Output arguments
varargout{1} = t_pre_spike_times;   % in ms
varargout{2} = tend;

end