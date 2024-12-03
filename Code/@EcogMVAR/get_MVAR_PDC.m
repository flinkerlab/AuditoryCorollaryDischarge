function [PDC, APDC] = get_MVAR_PDC(obj, A, varargin)
	% parse optional imports
	prs = inputParser;
	prs.KeepUnmatched = true;
	addParameter(prs, 'nf', 100);
	addParameter(prs, 'Q', []);
	addParameter(prs, 'mode', 1);
	%parse
	parse(prs, varargin{:});
	% set values
	nf = prs.Results.nf;
	Q = prs.Results.Q;
	mode = prs.Results.mode;
	Fs = obj.srate;

	% A is the size nd x nd x p x ns
	% nd: number of nodes
	% p : model order
	% ns: number of shift windows
	[nd,~,p,ns] = size(A);
	freqs = linspace(0,Fs/2,nf);
	APDC = zeros(nd,nd,nf,ns);
	Hf_temp = zeros(nd,nd,nf);
	if obj.flag_verbose; textprogressbar('Computing PDC:  '); end
	for s = 1:ns
		for f = 1:nf
			if mode==4
				if isempty(Q)
					[APDC(:,:,f,s),Hf_temp(:,:,f)] = EcogMVAR.pdc(A(:,:,:,s),...
																  eye(nd,nd),...
																  freqs(f),...
																  1/Fs, mode);
				else
					[APDC(:,:,f,s),Hf_temp(:,:,f)] = EcogMVAR.pdc(A(:,:,:,s),...
																  Q(:,:,s),...
																  freqs(f),...
																  1/Fs, mode);
				end
			else
				if mode==5
					if isempty(Q); 
						error('Q cannot be empty for mode=5!');
					else
						APDC(:,:,f,s) = EcogMVAR.pdc(A(:,:,:,s),...
													 Q(:,:,s),...
													 freqs(f),...
													 1/Fs, mode);

					end
				else
					APDC(:,:,f,s) = EcogMVAR.pdc(A(:,:,:,s),...
												 eye(nd,nd),...
												 freqs(f),...
												 1/Fs, mode);
				end
			end
		end	
		if mode==4
			DENTemp = (abs(APDC(:,:,:,s)).^2) .* (abs(Hf_temp).^2);
			APDC(:,:,:,s) = (abs(APDC(:,:,:,s)).*abs(Hf_temp))./(sum(DENTemp,'all')); 
		end
		if obj.flag_verbose; textprogressbar(s/ns*100); end
	end
	PDC = squeeze(sum(APDC,3));
	diaginds = [];
	for i=1:size(PDC,3)
        diaginds = [diaginds, [1:size(PDC,1)+1:size(PDC,1).^2]+(size(PDC,1).^2*(i-1))];
    end
    PDC(diaginds) = 0; % removing the self connections
	if obj.flag_verbose; textprogressbar(' done.'); end
end
