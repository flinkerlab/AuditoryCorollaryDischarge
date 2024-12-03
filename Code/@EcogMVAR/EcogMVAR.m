classdef EcogMVAR
	% class to fit and get AR model coefficents
	properties
		% model properties
		window_size = 20;
		window_shift = [];
		model_order = 3;
		srate = 200;
		flag_window_normalize = false;
		flag_verbose = true;
	end

	% methods for EcogMVAR class
	methods
		% Constructor function 
		function [obj] = EcogMVAR(varargin)
			% parse the input parameters
			p = inputParser;
			p.KeepUnmatched = true;
			addParameter(p, 'window_size', obj.window_size);
			addParameter(p, 'window_shift', floor(0.1*obj.window_size))
			addParameter(p, 'model_order', obj.model_order);
			addParameter(p, 'srate', obj.srate);
			addParameter(p, 'flag_window_normalize', obj.flag_window_normalize);
			addParameter(p, 'flag_verbose', obj.flag_verbose);
			% parse
			parse(p, varargin{:});
			% set parsed to obj
			obj.window_size = p.Results.window_size;
			obj.window_shift = p.Results.window_shift;
			obj.model_order = p.Results.model_order;
			obj.srate = p.Results.srate;
			obj.flag_window_normalize = p.Results.flag_window_normalize;
			obj.flag_verbose = p.Results.flag_verbose;
		end

		% fit MVAR function
		[A, Q, time_values, time_start_win] = FitMVAR(obj, TrialData, varargin);

		% wrap for pdc function
		[PDC, APDC] = get_MVAR_PDC(obj, A, varargin);
	end

	methods (Static)
		% the MEM fit function
		[X,Y,ind,aic,bic] = MEM(xt, order);
		
		% the fourier transform of A
		y = Af(A, f, dt);

		% Get pdc from a A for a window
		[y, Hout] = pdc(A, Q, f, dt, mode);

        % function to apply an AR model
        [x_hat, total_er, REV] = ARfilter(x, A, Q)

        % function to compute a val 
        [recon_sig, Val_ER, Val_REV] = ValidateMVAR(TrialData, A, Q, pt, s)

        % function to get max eig for the A comp
        [Eigvals] = get_MVAR_eigs(A, show_flag);
        
        % internal function to check MVAR fit on held-out trials
        [flag_check] = check_MVAR_coeffs(testTrials, A, Q, window_size, window_shift);

	end
end
