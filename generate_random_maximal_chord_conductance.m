function varargout = generate_random_maximal_chord_conductance(varargin)
% From Shouval 2004:
% The stochastic process for choosing G_NMDA has two phases:
%
% 1) We emulate the effect of failures of presynaptic release by randomly
% setting G_NMDA to zero, with a probability 0.5. It must be noted that, in
% our analysis, we do not differentiate between a release failure or a
% failure to bind postsynaptic receptors when there was a release.
%
% 2) Next, if there is a release, we chose the maximal chord conductance,
% G_NMDA randomly from a gamma distribution with a mean of 1/325, and a
% variance chosen differently for every \Delta t (spike timing difference)
% and the number of NMDA receptors, Z, according to linear fit to the
% theoretical derivation. The linear fits used in this paper are CV(\Delta
% t > 0) = 0.095 + 0.0045 \Delta t and CV(\Delta t <= 0) = 0.095 - 0.00067
% \Delta t. We also limit the maximal value of G_NMDA to be smaller than
% the mean NMDA receptor conductance times the number of NMDA receptors Z;
% this has a small effect on the results for small Z. Note that we do not
% change the conductance as a function of Z, only the relative variance.

% Input parameters
Delta_t = varargin{1};  % spike timing difference: t_post - t_pre
Z       = varargin{2};  % # of NMDA receptors

mu           = 1/325;     % mu M/(ms mV)
coeff_var    = 0.095 + 0.0045*max([0 Delta_t]) - 0.00067*min([0 Delta_t]); 
sigma_square = (mu*coeff_var)^2;
% --------------------------------------------------------------------
i = 3;
while i<=length(varargin),
    switch varargin{i},
        case 'sigma_square',    sigma_square = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% --------------------------------------------------------------------

% Model
if rand <= 0.5
    % Set G_NMDA to zero, with a probability 0.5
    G_NMDA = 0;
else
    % Choose it randomly from a gamma distribution, with a probability 0.5
    % under the condition that the maximal value of G_NMDA is smaller
    % than the mean NMDA receptor conductance times the number of NMDA
    % receptors Z
    option = 2;
    if option == 1
        while 1
            G_NMDA = gamrnd(mu^2/sigma_square, sigma_square/mu);
            if G_NMDA < mu*Z
                break;
            end
        end
    elseif option == 2
        while 1
            G_NMDA_vec = gamrnd(mu^2/sigma_square, sigma_square/mu, ...
                [2*Z 1]);
            G_NMDA_vec(G_NMDA_vec >= mu*Z) = [];
            if length(G_NMDA_vec) >= Z
                break;
            end
        end
        G_NMDA = mean(G_NMDA_vec(1:Z));
    end
end

% Output parameters
varargout{1} = (-1) * G_NMDA;   % Sign is set by Shouval 2002 

end