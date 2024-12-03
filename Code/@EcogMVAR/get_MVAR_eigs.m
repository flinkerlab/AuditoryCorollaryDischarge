function [Eigvals] = get_MVAR_eigs(A, show_flag)
	% gives the largest eigenvalue of the companion matrix for A
	% if A is nd x nd x p return a single number for 
	% if A is nd x nd x p x nw returns a vector of size nw 
	if nargin<2
		show_flag=false;
	end

	[nd,nd,p,nw] = size(A);
	I = speye((p-1)*nd);
	Z = sparse((p-1)*nd,nd);
	IZ = [I,Z];
	if nw==1
		A_comp = sparse(reshape(A,nd,nd*p));
		A_comp = [A_comp;IZ];
		Eigvals = abs(eigs(A_comp,1));
		return
	else
		Eigvals = zeros(nw,1);
		for i=1:nw
			A_comp = sparse(reshape(A(:,:,:,i),nd,p*nd));
			A_comp = [A_comp;IZ];
			Eigvals(i) = abs(eigs(A_comp,1));
		end
	end
	if show_flag
		figure;
		plot(1:nw,log10(Eigvals),'b*-', 'LineWidth',2, 'MarkerSize',7);
		ylim([min(log10(Eigvals))-0.1,0.1]);
		hold on;
		plot(1:nw,zeros(1,nw),'r--', 'LineWidth',2);
		grid on;
		xlabel('Window number'); ylabel('log10(max(Eig))');
	end
end

