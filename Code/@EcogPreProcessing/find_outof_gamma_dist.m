function [time_inds, thresh] = find_outof_gamma_dist(gdat, srate)
    % function to detect values out of dist from gamma
    % on analytic amplitude of high-gamma
    % fit one gamma to all the gdat after high-gamma filter
    gdat_hg = abs(EcogPreProcessing.hilbert_filter(gdat,srate,70,150));
    [p, ~] = gamfit(gdat_hg(:));
    values = 0:0.1:max(gdat_hg(:));
    gam_dist = gampdf(values, p(1), p(2));
    values = values(gam_dist>=1e-4);
    thresh = 30;%values(end);
    mask = zeros(size(gdat_hg));
    mask(gdat_hg>=thresh) = 1;
    time_inds = cell(size(gdat_hg,1),1);
    for i=1:size(gdat_hg,1)
        time_inds{i} = find(mask(i,:));
    end
end
