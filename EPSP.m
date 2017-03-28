function varargout = EPSP(varargin)

% Input arguments
t = varargin{1};                % Time (ms)
t_pre_spikes = varargin{2};     % Presynaptic spike times (ms)
norm = varargin{3};             % Chosen so that max EPSP is equal to s
which_parameters = varargin{4};

% Consider only the presynaptic spike that happened before t
t_pre_spikes = t_pre_spikes( t_pre_spikes < t );

% Fixed variables and functions
tau_1 = 50;
tau_2 = 5;

if which_parameters == 1
    % Shouval 2002: Figures 2, 3A and C, 5A and C
    s = 1;
elseif which_parameters == 2
    % Shouval 2002: Figures 3B and 5B
    s = 10;
elseif which_parameters == 3
    s = 20;    
end

% Output arguments
exp_terms = exp((t_pre_spikes-t)/tau_1) - exp((t_pre_spikes-t)/tau_2);
varargout{1} = s/norm * sum(exp_terms);

end
