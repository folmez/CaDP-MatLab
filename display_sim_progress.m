function varargout = display_sim_progress(varargin)
% [tDISP, tSIM] = display_sim_progress('initialize');

display_time_interval = 6; % seconds

switch varargin{1}
    case 'initialize'        
        % Initialize clocks
        tDISP = tic;     % Keep track display renewal
        tSIM  = tic;     % Keep track of the whole simulations
        
        % Output arguments
        varargout{1} = tDISP;
        varargout{2} = tSIM;
        
    otherwise
        % Input arguments
        i                = varargin{1};
        nr_time_steps    = varargin{2};
        tDISP            = varargin{3};
        tSIM             = varargin{4};
        display_progress = varargin{5};
        
        % Display simulation progress
        if display_progress
            if toc(tDISP) > display_time_interval
                % renew display clock
                tDISP = tic;
                % display progress
                fprintf('%2.2f%% is completed in %3.2f minutes...\n', ...
                    100*i/nr_time_steps, toc(tSIM)/60);
            end
        end
        
        % Output arguments
        varargout{1} = tDISP;
        
end

end