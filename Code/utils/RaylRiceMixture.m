classdef RaylRiceMixture
	% The class for Rayleigh-Rice Mixture Model Estimation
	properties
		% data
		Data
		% Rice properties
		RiceMixPortion = 0.5;
		RiceSigmaSqr = 1;
		RiceNu = 0;
		% Rayleigh properties
		RaylMixPortion = 0.5;
		RaylBSqr = 1;
		% Dist assignments
		RaylSamples = [];
		RiceSamples = [];
		% PDFs
		RicePDF = [];
		RaylPDF = [];
		MixturePDF = [];
		CalcPointPDF = [];
		% flags
		flag_verbose = true;
	end
	methods
		% Constructor Class;
		function [obj] = RaylRiceMixture(Data, varargin)
			% parse inputs
			p = inputParser;
			p.KeepUnmatched = true;
			addRequired(p, 'Data', @(d) (isnumeric(d)));
			addParameter(p, 'RiceMixPortion', obj.RiceMixPortion);
			addParameter(p, 'RiceSigmaSqr', obj.RiceSigmaSqr);
			addParameter(p, 'RiceNu', obj.RiceNu);
			addParameter(p, 'RaylMixPortion', obj.RaylMixPortion);
			addParameter(p, 'RaylBSqr', obj.RaylBSqr);
			addParameter(p, 'RaylSamples', obj.RaylSamples);
			addParameter(p, 'RiceSamples', obj.RiceSamples);
			addParameter(p, 'flag_verbose', obj.flag_verbose);
			addParameter(p, 'RicePDF', obj.RicePDF);
			addParameter(p, 'RaylPDF', obj.RaylPDF);
			addParameter(p, 'MixturePDF', obj.MixturePDF);
			addParameter(p, 'CalcPointPDF', obj.CalcPointPDF);
			%parse
			parse(p, Data, varargin{:});
			% collect the parameters
			obj.Data = p.Results.Data(:);
			% Rice properties
			obj.RiceMixPortion = p.Results.RiceMixPortion;
			obj.RiceSigmaSqr = p.Results.RiceSigmaSqr;
			obj.RiceNu = p.Results.RiceNu;
			% Rayleigh properties
			obj.RaylMixPortion = p.Results.RaylMixPortion;
			obj.RaylBSqr = p.Results.RaylBSqr;
			% Dist assignments
			obj.RaylSamples = p.Results.RaylSamples;
			obj.RiceSamples = p.Results.RiceSamples;
			% flags
			obj.flag_verbose = p.Results.flag_verbose;
			% PDF
			obj.RicePDF = p.Results.RicePDF;
			obj.RaylPDF = p.Results.RaylPDF;
			obj.MixturePDF = p.Results.MixturePDF;
			obj.CalcPointPDF = p.Results.CalcPointPDF;

		end

		% init function
		function [obj] = InitRaylRiceMixture(obj)
			% function to initilize the Mixture Model
			if any(obj.Data<0); error('Data cannot be negative!'); end
			x = obj.Data;
			[kind, kcnt] = kmeans(x, 2);
			if min(x(kind==1))<=min(x(kind==2))
				Raylind = 1; Riceind = 2;
			else
				Raylind = 2; Riceind = 1;
			end
			obj.RaylSamples = find(kind==Raylind);
			obj.RiceSamples = find(kind==Riceind);
			obj.RaylMixPortion = length(obj.RaylSamples)/length(x);
			obj.RiceMixPortion = length(obj.RiceSamples)/length(x);
			obj.RaylBSqr = raylfit(x(kind==Raylind), 0.05)^2;
			[nu,sig] = ricefit(x(kind==Riceind));
			obj.RiceNu = nu;
			obj.RiceSigmaSqr = sig^2;
			% calc PDFs
			obj.CalcPointPDF = linspace(0, max(x), 100);
			obj.RicePDF = obj.get_RicePDF(obj.CalcPointPDF,...
										  obj.RiceNu,...
										  obj.RiceSigmaSqr);
			obj.RaylPDF = obj.get_RaylPDF(obj.CalcPointPDF,...
										  obj.RaylBSqr);
			obj.MixturePDF = obj.RaylMixPortion .* obj.RaylPDF + ...
							 obj.RiceMixPortion .* obj.RicePDF;
		end

		function [obj] = FitRaylRiceMixture(obj, varargin)
			% function to fit the Mixture Model with EM algo
			prs = inputParser;
			prs.KeepUnmatched = true;
			addParameter(prs, 'tol', 1e-4);
			addParameter(prs, 'maxiter', 50);
			parse(prs, varargin{:});
			tol = prs.Results.tol;
			maxiter = prs.Results.maxiter;
			% init the mixture model
			obj = obj.InitRaylRiceMixture();
			% Take the values...
			alp_0 = obj.RaylMixPortion;
			BSqr_0 = obj.RaylBSqr;
			Nu_0 = obj.RiceNu;
			SigmaSqr_0 = obj.RiceSigmaSqr;
			X = obj.Data+1e-8;
			N = length(X);
			% start iteration
			if obj.flag_verbose; textprogressbar('Fitting RaylRice Mix Mdl: '); end
			for itr = 1:maxiter
				% find p(w|x,B,Nu,S)
				px1 = obj.get_RaylPDF(X, BSqr_0);
				px2 = obj.get_RicePDF(X, Nu_0, SigmaSqr_0);
				pw1 = alp_0 .* px1;
				pw2 = (1-alp_0) .* px2;
				pws = pw1+pw2;
				pw1 = pw1./pws;
				pw2 = pw2./pws;
				% update param
				alp_1 = (1/N) * sum(pw1,'all');
				BSqr_1 = sum(pw1.*(X.^2),'all') / (2*sum(pw1,'all'));
				I_1 = besseli(1, (X.*Nu_0)./(SigmaSqr_0));
				I_0 = besseli(0, (X.*Nu_0)./(SigmaSqr_0));
				Nu_1 = sum((pw2.*I_1.*X)./(I_0), 'all') ./ sum(pw2,'all');
				Temp = X.^2 + Nu_0.^2 - 2.*X.*Nu_0.*(I_1./I_0);
				SigmaSqr_1 = sum((pw2.*Temp), 'all') ./ (2*sum(pw2,'all'));
				% Check convergence
				if (abs(alp_0-alp_1)<=tol) && ...
				   (abs(BSqr_0-BSqr_1)<=tol) && ...
				   (abs(Nu_0-Nu_1)<=tol) && ...
				   (abs(SigmaSqr_0-SigmaSqr_1)<=tol);
						if obj.flag_verbose; textprogressbar('Converged!'); end
						break;
				else
					if obj.flag_verbose; textprogressbar(itr/maxiter*100); end
				end
				if itr==maxiter; 
					if obj.flag_verbose; textprogressbar('done!'); end
				end
				% move param for next iter
				alp_0 = alp_1;
				BSqr_0 = BSqr_1;
				Nu_0 = Nu_1;
				SigmaSqr_0 = SigmaSqr_1;
			end
			% set to obj
			obj.RaylMixPortion = alp_1;
			obj.RaylBSqr = BSqr_1;
			obj.RiceMixPortion = 1-alp_1;
			obj.RiceNu = Nu_1;
			obj.RiceSigmaSqr = SigmaSqr_1;
			% calc PDFs
			obj.CalcPointPDF = linspace(0, max(X), 100);
			obj.RicePDF = obj.get_RicePDF(obj.CalcPointPDF,...
										  obj.RiceNu,...
										  obj.RiceSigmaSqr);
			obj.RaylPDF = obj.get_RaylPDF(obj.CalcPointPDF,...
										  obj.RaylBSqr);
			obj.MixturePDF = obj.RaylMixPortion .* obj.RaylPDF + ...
							 obj.RiceMixPortion .* obj.RicePDF;
		end

		function [T] = get_RaylRiceMixThresh(obj, flag_Fit)
			if ~exist('flag_Fit', 'var'); flag_Fit = true; end
			if flag_Fit; obj = obj.FitRaylRiceMixture(); end
			% get values
			alp1 = obj.RaylMixPortion;
			BSqr = obj.RaylBSqr;
			alp2 = obj.RiceMixPortion;
			Nu = obj.RiceNu;
			SigmaSqr = obj.RiceSigmaSqr;
			x0 = (sqrt(BSqr)+Nu)/2;
			func = @(x) alp1*obj.get_RaylPDF(x,BSqr) -...
						alp2*obj.get_RicePDF(x,Nu,SigmaSqr);
			T = fzero(func, x0);
		end
		
	end

	% Static methods
	methods (Static)

		function [y] = get_RaylPDF(x, BSqr)
			% gives the raylpdf calculated on x
			y = (x./BSqr) .* (exp(-(x.^2)./(2*BSqr)));
			y(x<=0) = 0;
		end

		function [y] = get_RicePDF(x, Nu, SigmaSqr)
			% gives the raylpdf calculated on x
			y = (x./SigmaSqr) .*...
				(exp(-(x.^2+Nu.^2)./(2*SigmaSqr))) .*...
				besseli(0, (x.*Nu)./(SigmaSqr));
			y(x<=0) = 0;
		end
	end
end
