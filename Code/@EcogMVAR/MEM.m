function [ta, Sigma, ind_BIC, AIC, BIC] = MEM(xt, order)
    % function for fast estimation of model parameters
    % fits to model order selected
    % modified from supplemental material here:
    % https://www.eneuro.org/content/6/4/ENEURO.0472-18.2019.long
	[dim,ntr,nt] = size(xt);
	% The data needs to be dimension x number of time points x number of trials
	xt=permute(xt,[1 3 2]);
	a00 = 0;
	b00 = 0;
	p0 = 0;

	for r = 1:ntr
		%average over trials
		a00 = a00 + xt(:,(2:nt),r) * xt(:,(2:nt),r)';
		b00 = b00 + xt(:,1:(nt-1),r) * xt(:,1:(nt-1),r)';
		p0 = p0 + xt(:,:,r) * xt(:,:,r)';
	end

	Pn.p0 = chol( (1/(ntr*nt)) * p0,'lower');
	tr = Pn.p0;

	%average over trials, take cholesky decomposition and invert
	a.a00 = inv(chol( (1/ntr) * a00,'lower'));
	b.b00 = inv(chol( (1/ntr) * b00,'lower'));

	%initialize
	AIC = zeros(1,order);
	BIC = zeros(1,order);
	Sigma = cell(1,order);
	ta = cell(1,order);
	for n = 0:(order-1)
		a.(['a' num2str(n) num2str(n+1)]) = 0;
		b.(['b' num2str(n) num2str(n+1)]) = 0;

		%initialize
		rn_eps = 0;
		rn_r= 0;
		rn_epsr = 0;
		%average over trials
		for r = 1:ntr
			axt = 0;
			bxt = 0;
			for k = 0:n
				axt = axt + a.(['a' num2str(n) num2str(k)]) * (xt(:,(n-k+2):(nt-k),r));
				bxt = bxt + b.(['b' num2str(n) num2str(n-k)]) * (xt(:,(n+1-k):(nt-1-k),r));
			end
			%%%%%%%%%%
			rn_eps = rn_eps + axt * axt';
			rn_r = rn_r + bxt * bxt';
			rn_epsr = rn_epsr + axt * bxt';
		end

		%overwrite (only need inverses of rn_eps and rn_r)
		rn_eps = chol(rn_eps,'lower');
		rn_r = chol(rn_r,'lower');
		rho_n = (rn_eps \ rn_epsr) / rn_r';

		Pn.(['p' num2str(n+1)]) = (chol(eye(dim) - rho_n*rho_n','lower'));
		tr = tr * Pn.(['p' num2str(n+1)]);
		Sigma{n+1} = tr*tr';
		AIC(n+1) = log(det(Sigma{n+1}))*(nt*ntr) + 2*(dim^2)*(n+1);
		BIC(n+1) = log(det(Sigma{n+1}))*(nt*ntr) + (dim^2)*(n+1)*log(nt*ntr);


		Qn = (chol(eye(dim) - rho_n'*rho_n,'lower'));

		tmpa = [];
		for k = 0:(n+1)

			nk = [num2str(n) num2str(k)];
			nk_1 = [num2str(n) num2str(n-k+1)];
			a.(['a' num2str(n+1) num2str(k)]) = Pn.(['p' num2str(n+1)]) \ (a.(['a' nk]) - rho_n*b.(['b' nk_1]));
			b.(['b' num2str(n+1) num2str(k)]) = Qn \ (b.(['b' nk]) - rho_n'*a.(['a' nk_1]));

			if(k>0)
				tmpa = [tmpa, (a.(['a' num2str(n+1) '0'])) \ a.(['a' num2str(n+1) num2str(k)])];
			end
		end
		ta{n+1}=tmpa;
	end
	[minVal,ind_AIC] = min(AIC);
	[minVal,ind_BIC] = min(BIC);
    % plot(AIC+BIC); hold on;

	%parameter specified in papers (should be < 0.1)
	paramOrd = dim * (ind_BIC+1) / (nt*ntr);
	if(paramOrd > 0.1)
		fprintf('warning ParamOrd:%2.2f, indBIC:%2d \n',paramOrd,ind_BIC);
	end
end
