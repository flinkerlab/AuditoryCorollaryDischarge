function [Acnct, cnct] = get_flatten_connectivity(A);
    % reshape the connectivty matrix into cnct x time
    % A -- matrix of size:
    % -- coeff: node x node x model_order x time
    % -- pdc:   node x node x time

    if length(size(A)) == 4
        A = squeeze(sum(abs(A),3));
    end
    if (size(A,1)~=size(A,2));
        error('matrix A needs to be square in 1st and 2nd dims!');
    end
    n = size(A,1);
    t = size(A,3);
    Acnct = reshape(A,n.^2,t);
    Acnct(1:n+1:n^2,:) = [];

    clist = 1:n^2;
    clist(1:n+1:n^2) = [];
    [cnct.from, cnct.to] = ind2sub([n,n], clist);
    cnct.cnct_list = clist;
end
