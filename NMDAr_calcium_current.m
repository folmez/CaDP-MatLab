function varargout = NMDAr_calcium_current(varargin)

% Input arguments
t = varargin{1};            % Time
t_pre_spikes = varargin{2}; % Presynaptic spike time (in ms)
V   = varargin{3};          % Postysnaptic membrane potential (in mV)
closed_NMDAr_frac_before_spikes = varargin{4};
version = 'deterministic';      % 'deterministic' (Shouval 2002)
                                % 'deterministic-with-stochastic-parameters'
                                % 'stochastic'    (Shouval 2004)
tau_f  = 50;                % Fast component decay rate (in ms)
tau_s  = 200;               % Slow component decay rate (in ms)
% -------------------------------------------------------------------------
i = 5;
while i<=length(varargin),
    switch varargin{i},
        case 'G_NMDA_for_each_spike'
            G_NMDA_for_each_spike = varargin{i+1};
        case 'version',         version = varargin{i+1};
        case 'NMDAr_I_f',           I_f = varargin{i+1};
        case 'NMDAr_tau_f',       tau_f = varargin{i+1};
        case 'NMDAr_tau_s',       tau_s = varargin{i+1};
        case 'spt',                 spt = varargin{i+1};
        case 'Z',                     Z = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% -------------------------------------------------------------------------
% Consider only the presynaptic spike that happened before t
t_pre_spikes = t_pre_spikes( t_pre_spikes < t );
t_i = t-t_pre_spikes;

% Fixed variables and functions
if strcmp(version, 'deterministic')
    P0     = 0.5;       % Fraction of NMDArs that move from closed to open
    G_NMDA = -1/500;    % Peak NMDAr conductance (in muM/(ms*mV))
else
    I_f    = 0.75;      % Relative magnitude of the fast component
    tau_f  = 50;
    tau_s  = 150;
    if strcmp(version, 'deterministic-with-stochastic-parameters')
        P0     = 0.5;
        G_NMDA = -1/325; % this is + in Shouval 2004, I think it is typo
    elseif strcmp(version, 'stochastic')
        P0     = 0.5;
        G_NMDA = G_NMDA_for_each_spike;
    end
end
V_r    = 130;       % Reversal potential for calcium (in mV)
I_s    = 1-I_f;     % Relative magnitude of the slow component

Mg = 1;                                   % Magnesium (Mg) concentration
B  = @(V) 1/(1+exp(-0.062*V)*(Mg/3.57));  % Effect of Mg block    

% Output arguments
if ~isempty(t_i) && t_i(end) > 0
    NMDAr_frac_vec = closed_NMDAr_frac_before_spikes .* P0 .* ...
        (I_f*exp(-t_i/tau_f) + I_s*exp(-t_i/tau_s));
    
    total_NMDAr_frac = sum(NMDAr_frac_vec);
    varargout{1} = sum( (NMDAr_frac_vec.*G_NMDA) ) * B(V) * (V-V_r);
    varargout{2} = total_NMDAr_frac;
else
    varargout{1} = 0;
    varargout{2} = 0;
end

end

% Previous calculation method
% t_i = t_i(end);
% varargout{1} = P0 * G_NMDA * ...
%     ( I_f*exp(-t_i/tau_f) + I_s*exp(-t_i/tau_s) ) * B(V) * (V-V_r);
% varargout{2} = P0 * ( I_f*exp(-t_i/tau_f) + I_s*exp(-t_i/tau_s) );