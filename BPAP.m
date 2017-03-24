function varargout = BPAP(varargin)
% BPAP is the back-propagating action potential. 100 is the maximal value
% in mV that the BPAP will depolarize the cell to above resting potential.
% In Shouval 2002, simulations except in Figs. 2A and 4A use the same set
% of parameters: tau_f = 3 ms , tau_s = 25 ms. In Fig. 4A we used tau_f =
% 15 ms and tau_s = 50 ms . A 2-ms delay is assumed in the arrival time of
% the BPAP to the spine.

% Input arguments
t = varargin{1};
t_post_spike = varargin{2};
BPAP_delay = varargin{3};
BPAP_type = varargin{4};

% BPAP arrives late, the spine should receive the BPAP from before
t = t-BPAP_delay;

% Consider only the postsynaptic spike that happened after t
t_post_spike = t_post_spike(t_post_spike<t);
if isempty(t_post_spike)
    % no postsynaptic spike occurred yet
    t_post_spike = -inf;
else
    % take the latest postsynaptic spike
    t_post_spike = t_post_spike(end);
end

% Fixed variables and functions
I_f   = 0.75;   % Relative magnitude of the fast component
I_s   = 1-I_f;  % Relative magnitude of the slow component
if strcmp(BPAP_type, 'BPAP + ADP')
    % ADP: After-depolarizing potential, i.e. slow component
    tau_f = 3;      % Fast component decay rate (in ms)
    tau_s = 25;     % Slow component decay rate (in ms)
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
if isempty(t_post_spike)
    varargout{1} = 0;
else
    varargout{1} = 100 * ( I_f * exp(-(t-t_post_spike)/tau_f) + ...
        I_s * exp(-(t-t_post_spike)/tau_s));
end

end
