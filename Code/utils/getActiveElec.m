function [ActiveInd] = getActiveElec(X, varargin)
    p = inputParser;
    addParameter(p, 'mode', 'smoothRR');
    addParameter(p, 'thresh', 1);
    addParameter(p, 'verbose', false);
    parse(p, varargin{:});
    options = p.Results;
    
    % different active elec detection algo
    if options.verbose
        fprintf('detecting active elec...\n');
    end
    if strcmp(options.mode, 'smooth')
        xm = squeeze(mean(x,2));
        xm_smooth = zeros(size(xm));
        for i = 1:size(xm,1)
			xm_smooth(i,:) = WTDN(xm(i,:), 'db8', 0.5);
		end
		sd_all = std(xm_smooth, [], 2);
		[rice_mean,rice_sigma] = ricefit(sd_all);
		ActiveInd = find(sd_all>=options.thresh*rice_sigma);
	elseif strcmp(options.mode, 'smoothRR')
		xm = squeeze(mean(x,2));
		xm_smooth = zeros(size(xm));
		for i = 1:size(xm,1)
			xm_smooth(i,:) = WTDN(xm(i,:), 'db8', 0.5);
		end
		sd_all = std(xm_smooth, [], 2);
		RR = RaylRiceMixture(sd_all);
		thresh = RR.get_RaylRiceMixThresh();
		ActiveInd = find(sd_all>=thresh);
	elseif strcmp(options.mode, 'Fixed')
		xm = squeeze(mean(x,2));
		xm_smooth = zeros(size(xm));
		for i = 1:size(xm,1)
			xm_smooth(i,:) = WTDN(xm(i,:), 'db8', 0.5);
		end
		sd_all = std(xm_smooth, [], 2);
		ActiveInd = find(sd_all>=options.thresh);
	else
		error('mode not implemented!')
	end
end
