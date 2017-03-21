function varargout = NMDAr_calcium_current(varargin)

% Input arguments
t = varargin{1};            % Time
t_pre_spike = varargin{2};  % Presynaptic spike time (in ms)
V   = varargin{3};          % Postysnaptic membrane potential (in mV)

t_i = t-t_pre_spike;

% Fixed variables and functions
P0     = 0.5;       % Fraction of NMDArs that move from closed to open
G_NMDA = -1/500;    % Peak NMDAr conductance (in muM/(ms*mV))
V_r    = 130;       % Reversal potential for calcium (in mV)
I_f    = 0.5;
I_s    = 1-I_f;     % Relative magnitude of the slow component
tau_f  = 50;        % Fast component decay rate (in ms)
tau_s  = 200;       % Slow component decay rate (in ms)

Mg = 1;                                     % Magnesium (Mg) concentration
B  = @(V) 1/(1+exp(-0.062*V)*(Mg/3.57));  % Effect of Mg block    

% Output arguments
if t_i>0
    varargout{1} = P0 * G_NMDA * ...
        ( I_f*exp(-t_i/tau_f) + I_s*exp(-t_i/tau_s) ) * B(V) * (V-V_r);
else
    varargout{1} = 0;
end
if length(varargout)>1
    varargout{2} = P0 * heaviside(t_i) * ...
        ( I_f*exp(-t_i/tau_f) + I_s*exp(-t_i/tau_s) );
end

end