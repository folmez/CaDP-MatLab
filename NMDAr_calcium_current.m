function varargout = NMDAr_calcium_current(varargin)

% Input arguments
t = varargin{1};            % Time
t_pre_spikes = varargin{2}; % Presynaptic spike time (in ms)
V   = varargin{3};          % Postysnaptic membrane potential (in mV)
I_f = varargin{4};

which_calculation_method = 2;

t_i = t-t_pre_spikes;

% Fixed variables and functions
P0     = 0.5;       % Fraction of NMDArs that move from closed to open
G_NMDA = -1/500;    % Peak NMDAr conductance (in muM/(ms*mV))
V_r    = 130;       % Reversal potential for calcium (in mV)
I_s    = 1-I_f;     % Relative magnitude of the slow component
tau_f  = 50;        % Fast component decay rate (in ms)
tau_s  = 200;       % Slow component decay rate (in ms)

Mg = 1;                                   % Magnesium (Mg) concentration
B  = @(V) 1/(1+exp(-0.062*V)*(Mg/3.57));  % Effect of Mg block    

% Output arguments
if ~isempty(t_i) && t_i(end) > 0
    if which_calculation_method == 1
        t_i = t_i(end);
        varargout{1} = P0 * G_NMDA * ...
            ( I_f*exp(-t_i/tau_f) + I_s*exp(-t_i/tau_s) ) * B(V) * (V-V_r);
        varargout{2} = P0 * ( I_f*exp(-t_i/tau_f) + I_s*exp(-t_i/tau_s) );
    elseif which_calculation_method == 2
        closed_NMDAr_frac_before_spikes = varargin{5};
        NMDAr_frac = sum( closed_NMDAr_frac_before_spikes .* P0 .* ...
            (I_f*exp(-t_i/tau_f) + I_s*exp(-t_i/tau_s)) );
        varargout{1} = NMDAr_frac * G_NMDA * B(V) * (V-V_r);
        varargout{2} = NMDAr_frac;
    end
else
    varargout{1} = 0;
    varargout{2} = 0;
end

end