function [recon_sig, Val_ER, Val_REV] = ValidateMVAR(TrialData, A, Q, pt, s)
	% Validation error of the MVAR model the trail data
	% TrialData: nd x tr x nb
	% pt: quasi-stationary window size
	% s: window shift 

	if nargin<5; s = floor(0.1*pt); end
	
	[nd,tr,nb] = size(TrialData);
	N_windows = floor((nb-pt)/s); % number of quasi-stationary windows
	Val_ER = zeros(1,N_windows);
	Val_REV = zeros(1,N_windows);
    recon_sig = nan*ones(size(TrialData));
	for n = 0:N_windows
		[temp,Val_ER(n+1),Val_REV(n+1)] = EcogMVAR.ARfilter(TrialData(:,:,s*n+1:s*n+pt),...
												A(:,:,:,n+1), Q(:,:,n+1));
        recon_sig(:,:,s*n+1:s*n+pt) = temp;
	end
end
