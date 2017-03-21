function varargout = calculate_spike_times(varargin)

% Input arguments
pre_spike_freq = varargin{1};    % Presynaptic freq (in Hz)
nr_pre_spikes  = varargin{2};    % # of presynaptic pulses
t0             = varargin{3}(1); % in ms
tend           = varargin{3}(2); % in ms

t_unit = 1000;                   % s to ms conversion factor

% Check inputs
if nr_pre_spikes/((tend-t0)/t_unit) > pre_spike_freq
    error('FO: Time interval too short for presynaptic stimulation.');
end
    
% Model
inter_spike_int = t_unit/pre_spike_freq;                       % in ms
t_pre_spike_times = t0 + (0:nr_pre_spikes-1)'*inter_spike_int; % in ms

% Output arguments
varargout{1} = t_pre_spike_times;   % in ms

end