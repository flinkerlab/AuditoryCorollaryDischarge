function [Wa, W1, W2] = get_clustering_weights(U,Ulb);
    % form the Adjacency matrix per component
    W1 = zeros(size(U,[2,3]));
    W2 = zeros(size(U,[2,3]));
    Wa = zeros(size(U,[2,3]));
    n_cnct_comp = zeros(size(U,3),1);
    n_cnct_totl = size(U,1)^2-size(U,1);

    for k = 1:size(U,3);
        n_cnct_comp(k) = length(find(Ulb==k));
        W1(:,k) = sum(squeeze(U(:,:,k)),1)./n_cnct_comp(k);
        W2(:,k) = sum(squeeze(U(:,:,k)),2)./n_cnct_comp(k);
    end
    Wa = W2-W1;
end
