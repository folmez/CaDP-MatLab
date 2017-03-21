function varargout = update_Ca(varargin)

% Input parameters
Ca           = varargin{1};     % Previous calcium level
NMDA_cal_cur = varargin{2};     % NMDA calcium current
dt           = varargin{3};     % Time-step (ms)

% Fixed variables and functions
tau_Ca = 50;    % in ms

% Output parameters
dCa_over_dt = NMDA_cal_cur - (1/tau_Ca)*Ca;
varargout{1} = Ca + dCa_over_dt * dt;

end