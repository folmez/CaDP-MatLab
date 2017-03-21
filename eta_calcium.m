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

    % Fixed variables and functions
    P1 = 0.1;       % in seconds
    P2 = P1*1e-4;   % unitless (in Shouval 2002, this parameter has a typo)
    P3 = 3;         % in seconds
    P4 = 1;         % in seconds

    % Output arguments    
    tau = P1 ./ (P2+Ca.^P3) + P4;
    varargout{1} = 1./tau;
end

end

function test_eta_calcium
Ca_vector = linspace(0, 1, 101)';
eta_calcium_vals = eta_calcium(Ca_vector);

figure,
plot(Ca_vector, eta_calcium_vals);
xlabel('Ca^2^+ (\muM)');
ylabel('\eta (sec^-^1)');
ylim([-0.25 1.25]);
title('\eta - Ca-dependent learning rate', 'FontSize', 15);
end