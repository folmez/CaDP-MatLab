function varargout = Omega_calcium(varargin)

if strcmp(varargin{1}, 'test')
    test_Omega_calcium();
else
    % Input arguments
    Ca = varargin{1};   % Calcium level
    
    % Fixed variables and functions
    sig = @(x, beta) exp(beta*x)./(1+exp(beta*x));
    alpha1 = 0.35;
    alpha2 = 0.55;
    beta1  = 80;
    beta2  = 80;
    
    % Output arguments
    varargout{1} = 0.25 + sig(Ca-alpha2, beta2) - ...
        0.25*sig(Ca-alpha1, beta1);
end

end

function test_Omega_calcium
Ca_vector = linspace(0, 1, 101)';
Omega_calcium_vals = Omega_calcium(Ca_vector);

figure,
plot(Ca_vector, Omega_calcium_vals);
xlabel('Ca^2^+ (\muM)');
ylabel('\Omega');
ylim([-0.25 1.25]);
hold on;
plot([0 1], [0.25 0.25], 'k--');
end