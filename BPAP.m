function varargout = BPAP(varargin)
% BPAP is the back-propagating action potential. 100 is the maximal value
% in mV that the BPAP will depolarize the cell to above resting potential.
% In Shouval 2002, simulations except in Figs. 2A and 4A use the same set
% of parameters: tau_f = 3 ms , tau_s = 25 ms. In Fig. 4A we used tau_f =
% 15 ms and tau_s = 50 ms . A 2-ms delay is assumed in the arrival time of
% the BPAP to the spine.

% Input arguments
t = varargin{1};
t_post_spikes = varargin{2};
BPAP_delay = varargin{3};
BPAP_type = varargin{4};
BPAP_tau_s = 25;                % BPAP slow decay rate
% --------------------------------------------------------------------
i = 5;
while i<=length(varargin),
    switch varargin{i},
        case 'BPAP_tau_s',            BPAP_tau_s = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% --------------------------------------------------------------------
output_option = 1;

% BPAP arrives late, the spine should receive the BPAP from before
t = t-BPAP_delay;

% Consider only the postsynaptic spike that happened after t
t_post_spikes = t_post_spikes( t_post_spikes < t );
if isempty(t_post_spikes)
    % no postsynaptic spike occurred yet
    last_t_post_spike = -inf;
else
    % take the latest postsynaptic spike
    last_t_post_spike = t_post_spikes(end);
end

% Fixed variables and functions
I_f   = 0.75;   % Relative magnitude of the fast component
I_s   = 1-I_f;  % Relative magnitude of the slow component
if strcmp(BPAP_type, 'BPAP + ADP')
    % ADP: After-depolarizing potential, i.e. slow component
    tau_f = 3;              % Fast component decay rate (in ms)    
    tau_s = BPAP_tau_s;     % Slow component decay rate (in ms)   
elseif strcmp(BPAP_type, 'not set yet')
    tau_f = 15;
    tau_s = 50;
elseif strcmp(BPAP_type, 'Narrow BPAP')
    % Only the fast component
    tau_f = 3;
    I_f   = 1;
    tau_s = 1;
    I_s   = 0;
end

% Output arguments
if output_option == 1
    varargout{1} = 100 * ( I_f * exp(-(t-last_t_post_spike)/tau_f) + ...
        I_s * exp(-(t-last_t_post_spike)/tau_s));
elseif output_option == 2
    varargout{1} = 100 * ( I_f * sum(exp(-(t-t_post_spikes)/tau_f)) + ...
        I_s * sum(exp(-(t-t_post_spikes)/tau_s)));
end

end
