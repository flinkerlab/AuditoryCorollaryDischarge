function [A, Q, time_values, time_start_win] = FitMVAR(obj, TrialData, varargin)
	% Fit the MVAR model to the trail data
	% TrialData: nd x tr x nb
	% np: model order
	% pt: quasi-stationary window size
	% s: window shift 
    argin = inputParser;
    addParameter(argin, 'time_vector', []);
    parse(argin, varargin{:});
    time_vector = argin.Results.time_vector;

	% get sizes
	[nd,tr,nb] = size(TrialData);
	% number of quasi-stationary windows
	N_windows = floor((nb-obj.window_size)/obj.window_shift); 
	A = zeros(nd,nd,obj.model_order,N_windows+1);
	Q = zeros(nd,nd,N_windows+1);
	pt = obj.window_size;
	s  = obj.window_shift;
    if ~isempty(time_vector);
        time_values = zeros(N_windows+1, pt);
        time_start_win = zeros(N_windows+1,1);
    end
	if obj.flag_verbose; textprogressbar('Fitting MVAR model: '); end
	for n = 0:N_windows
		if obj.flag_window_normalize
			xx = raw_window_preproc(TrialData(:,:,s*n+1:s*n+pt));
		else
			xx = TrialData(:,:,s*n+1:s*n+pt);
		end
		[X,Y,ind,aic,bic] = obj.MEM(xx, obj.model_order);
		A(:,:,:,n+1) = reshape(X{obj.model_order},nd,nd,obj.model_order);
		Q(:,:,n+1) = Y{obj.model_order};
        if ~isempty(time_vector)
            time_values(n+1,:) = time_vector(s*n+1:s*n+pt);
            time_start_win(n+1) = time_vector(s*n+1);
        end
		if obj.flag_verbose; textprogressbar(n/N_windows*100); end
	end
	if obj.flag_verbose; textprogressbar('done.'); end

	% time for each window
    if isempty(time_vector)
        time_values = {};
        time_start_win = [];
        shift_window_time = s * (1/obj.srate) * 1000;
        full_window_time = pt * (1/obj.srate) * 1000;
        for i=1:N_windows+1
            time_start_win(end+1) = (i-1)*shift_window_time;
            time_values{end+1} = [num2str((i-1)*shift_window_time), '-',...
            num2str((i-1)*shift_window_time+full_window_time)];
        end
    end
end
