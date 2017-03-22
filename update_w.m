function varargout = update_w(varargin)

% Input parameters
w  = varargin{1};     % Previous weight
Ca = varargin{2};     % Calcium level
dt = varargin{3};     % Time-step (ms)
lambda = 1;

% Optional inputs
i = 4;
while i<=length(varargin),
    switch varargin{i},
        case 'lambda',  lambda= varargin{i+1};
        otherwise,      display(varargin{i}); error('Unexpected inputs!!!');
    end
    i = i+2;
end

% Fixed variables and functions
t_unit = 1000;        % sec to ms conversion factor (to be used for eta) 

% Output parameters
dw_over_dt = eta_calcium(Ca)/t_unit * (Omega_calcium(Ca) - lambda * w);
varargout{1} = w + dw_over_dt * dt;

end