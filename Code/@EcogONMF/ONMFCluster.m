function [V, Ulb, U] = ONMFCluster(obj, A)
	% get connection matrix
	[Aflat, cnct] = obj.get_flatten_connectivity(A);

	% Perform Non-negative Matrix Factorization
	% mode == "onmf" --> orthogonal nnmf
	% mode == "kmeans" --> kmeans
	if strcmp(obj.mode, "onmf") % orth nnmf multiple runs
		Rsq = realmax;
		V = []; U_lbls=[];
		if obj.flag_verbose
			textprogressbar('Clustering with ONMF: ');
		end
		% iterate the clustering and choose the min relError
		for iter = 1:obj.maxiter
			[U_iter, V_iter, relError] = obj.FitONMF(Aflat,...
													 obj.num_clusters,...
													 obj.Options.ONMF_inner_iters,...
													 0);
			if (relError)<=Rsq
				Rsq = (relError); V = V_iter; U_lbls = U_iter;
			end
			if obj.flag_verbose
				textprogressbar(iter/obj.maxiter*100);
			end
		end
        % form the U matrix
        U = zeros([size(Aflat,1), obj.num_clusters]);
        for k_ind = 1:obj.num_clusters
            idx_k = (U_lbls==k_ind);
            Mi = Aflat(idx_k,:);
            U(idx_k,k_ind) = Mi*V(k_ind,:)';
        end
		% end of clustering
		if obj.flag_verbose
			textprogressbar('done.');
		end
		% if requested, order the clusters...
		% order based on temporal peak
		if strcmp(obj.Options.flag_order_centroids, 'time')
			if obj.flag_verbose
				fprintf('ordering clusters based on time ...');
			end
			[MaxVals, MaxInds] = max(V, [], 2);
			[MaxIndsSorted, sind] = sort(MaxInds);
			V_sorted = zeros(size(V));
			U_sorted = zeros(size(U));
			U_lbls_sorted = zeros(size(U_lbls));
			for j = 1:length(sind)
				U_lbls_sorted(U_lbls==sind(j)) = j;
				V_sorted(j,:) = V(sind(j),:);
				U_sorted(:,j) = U(:,sind(j));
			end
			V = V_sorted; U= U_sorted;
			U_lbls = U_lbls_sorted;
			if obj.flag_verbose
				fprintf(' done.\n');
			end
		end
		% reshape the clusters and weights to connection matrix
		Ulb = zeros(size(A,1));
        Ulb(cnct.cnct_list) = U_lbls;
        % reshape U to a cnct shape matrix
        U_mat = zeros([size(A,1), size(A,1), obj.num_clusters]);
        for i=0:obj.num_clusters-1
            U_mat(cnct.cnct_list+(i*size(A,1)^2)) = U(:,i+1);
        end
        U = U_mat;
	else
		error('Not implemented mdoe');
	end
end

