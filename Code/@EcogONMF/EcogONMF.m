classdef EcogONMF
	% properties for the EcogONMF
	properties
		mode = "onmf";
		maxiter = 20;
		num_clusters = 5;
		Options = struct();
		flag_verbose = true;
	end

	% methods
	methods
		% constructor function
		function [obj] = EcogONMF(varargin)
			% constructor function
			p = inputParser;
			p.KeepUnmatched = true;
			addParameter(p, 'mode', "onmf");
			addParameter(p, 'maxiter', obj.maxiter);
			addParameter(p, 'num_clusters', obj.num_clusters);
			addParameter(p, 'ONMF_inner_iters', 500);
			addParameter(p, 'flag_order_centroids', 'time');
			addParameter(p, 'flag_Diagnose', true);
			addParameter(p, 'flag_verbose', obj.flag_verbose);
			%parse
			parse(p, varargin{:});
			% set parameters
			obj.mode = p.Results.mode;
			obj.maxiter = p.Results.maxiter;
			obj.num_clusters = p.Results.num_clusters;
			obj.Options.ONMF_inner_iters = p.Results.ONMF_inner_iters;
			obj.Options.flag_order_centroids = p.Results.flag_order_centroids;
			obj.Options.flag_Diagnose = p.Results.flag_Diagnose;
			obj.flag_verbose = p.Results.flag_verbose;
		end

		% the clustering function
		[V, Ulb, U] = ONMFCluster(obj, A);

	end

	% static methods
	methods (Static)
		[labels,...
		 centroids,...
		 relError,...
		 iter,...
		 U_binarized] = FitONMF(M,...
								initClusts,...
								MaxIters,...
								initMethod);

        % function to flatten A or PDC
        [Acnct, cnct] = get_flatten_connectivity(A);

        % function to get weights for plotting
        [Wa, W1, W2] = get_clustering_weights(U,Ulb);

        % function to visualize the results
        visual_clustering(V, Ulb, U, coords, varargin);

        visual_clustering_suppression(V, Ulb, U, coords, varargin);
	
	end
end
