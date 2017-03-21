function varargout = update_w(varargin)

% Input parameters
w  = varargin{1};     % Previous weight
Ca = varargin{2};     % Calcium level
dt = varargin{3};     % Time-step (ms)

% Fixed variables and functions
t_unit = 1000;        % sec to ms conversion factor (to be used for eta) 

% Output parameters
dw_over_dt = eta_calcium(Ca)/t_unit * (Omega_calcium(Ca) - w);
varargout{1} = w + dw_over_dt * dt;

end