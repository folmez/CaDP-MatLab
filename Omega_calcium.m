function varargout = Omega_calcium(varargin)

if strcmp(varargin{1}, 'test')    
    test_Omega_calcium();
else
    % Input arguments
    Ca      = varargin{1};          % Calcium level
    version = 'deterministic';      % 'deterministic' (Shouval 2002)
                                    % 'stochastic'    (Shouval 2004)
    % --------------------------------------------------------------------
    i = 2;
    while i<=length(varargin),
        switch varargin{i},
            case 'version'
                version = varargin{i+1};
            otherwise,
                display(varargin{i});
                error('Unexpected inputs!!!');
        end
        i = i+2;
    end
    % --------------------------------------------------------------------
    % Fixed variables and functions
    sig = @(x, beta) exp(beta*x)./(1+exp(beta*x));    
    A = [0.25 1.00 -0.25];  % A(3) has a typo in Shouval 2004
    if strcmp(version, 'deterministic')
        alpha1 = 0.35;      alpha2 = 0.55;
        beta1  = 80;        beta2  = 80;
    elseif strcmp(version, 'stochastic')
        alpha1 = 0.40;      alpha2 = 0.65;
        beta1  = 30;        beta2  = 30;
    end
    
    % Output arguments
    varargout{1} = A(1) + A(2) * sig(Ca-alpha2, beta2) + ...
        A(3) * sig(Ca-alpha1, beta1);
end

end

function test_Omega_calcium
Ca_vector = linspace(0, 1, 101)';
Omega_ca_det_vals = Omega_calcium(Ca_vector, 'version', 'deterministic');
Omega_ca_sto_vals = Omega_calcium(Ca_vector, 'version', 'stochastic');

figure,
plot(Ca_vector, Omega_ca_det_vals);
hold on;
plot(Ca_vector, Omega_ca_sto_vals);
xlabel('Ca^2^+ (\muM)');
ylabel('\Omega');
% ylim([-0.25 1.25]);
ylim([0 1.05]);
hold on;
plot([0 1], [0.25 0.25], 'k--');
legend('Deterministic (Shouval 2002)', 'Stochastic (Shouval 2004)', ...
    'Location', 'Best');
end