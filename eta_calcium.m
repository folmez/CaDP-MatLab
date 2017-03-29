function varargout = eta_calcium(varargin)
% ETA_CALCIUM is the calcium-dependent learning rate. The parameters P1-P4
% are set so that when [Ca] ~ 0, tau ~ 3h. Implicitly throughout the paper
% (Shouval 2002), it is assumed that resting calcium levels are zero. This
% implies that calcium is measured with respect to its resting value, which
% are believed to be 50–100 nM.

% IMPORTANT NOTE: eta is 1/sec, needs to be converted to ms

if strcmp(varargin{1}, 'test')
    test_eta_calcium();
else
    % Input arguments
    Ca = varargin{1};   % Calcium level
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
    if strcmp(version, 'deterministic')
        P1  = 0.1;       % in seconds
        P2  = P1*1e-4;   % unitless (in Shouval 2002, there is a typo)
        P3  = 3;         % in seconds
        P4  = 1;         % in seconds        
        tau = P1 ./ (P2+Ca.^P3) + P4;
    elseif strcmp(version, 'stochastic')
        P1  = 1;
        P2  = 0.6;
        P3  = 3;
        P4  = 0.00001;
        tau = 1./( (P1*(Ca+P4).^P3) ./ ((Ca+P4).^P3+P2^P3) );
    end
    
    % Output arguments    
    varargout{1} = 1./tau;
end

end

function test_eta_calcium
Ca_vector = linspace(0, 1, 101)';
eta_calcium_det_vals = eta_calcium(Ca_vector, 'version', 'deterministic');
eta_calcium_sto_vals = eta_calcium(Ca_vector, 'version', 'stochastic');

figure,
plot(Ca_vector, eta_calcium_det_vals);
hold on;
plot(Ca_vector, eta_calcium_sto_vals);
xlabel('Ca^2^+ (\muM)');
ylabel('\eta (sec^-^1)');
ylim([-0.25 1.25]);
title('\eta - Ca-dependent learning rate', 'FontSize', 15);
legend('Deterministic (Shouval 2002)', 'Stochastic (Shouval 2004)', ...
    'Location', 'Best');
end