function [filtered_signal] = meyer_multiband(gdat, Fs, fl, fh, n_bands)
    % create a uniform frequency grid
    n_samples = size(gdat,2);
    L_fft = 2^nextpow2(n_samples);
    freq = linspace(0, Fs/2, L_fft/2+1);

    % compute band-width
    alpha = 1.5;
    transit_width = 5;
    total_band_width = fh-fl;
    % df = total_band_width/(2^(n_bands+1)-1); % smallest band-width
    df = total_band_width*(alpha-1)/(alpha^(n_bands+1)-1);

    % compute the middle points of the bands [fm-2^i*df, fm+2^(i+1)*df]
    % f_middle = (2.^(1:n_bands)-1)*df + fl;
    f_middle = (alpha.^(1:n_bands)-1)/(alpha-1)*df + fl;
    filter_bank = zeros(n_bands+2, L_fft/2+1);

    % 0th band is [fl-df, fl+df] symmetric --> not in the f_middle loop
    % i-th band is [f_middle(i)-2^(i-1)*df, f_middle(i)+2^i*df]
    % last band is [f_middle(end), fh+df] --> sine up with band fh-f_middle(end), cos down df

    lb = fl-transit_width; ub = fl+df;
    ind_lb = find((freq>=lb)&(freq<=fl));
    filter_bank(1,ind_lb) = sin((pi/2)*meyeraux((freq(ind_lb)-lb)/transit_width));
    ind_ub = find((freq>=fl)&(freq<=ub));
    filter_bank(1,ind_ub) = cos((pi/2)*meyeraux((freq(ind_ub)-fl)/df));

    % loop over bands
    for i = 1:n_bands
        lb = f_middle(i) - alpha^(i-1)*df; ub = f_middle(i) + 2^i*df;
        ind_lb = find((freq>=lb)&(freq<=f_middle(i)));
        filter_bank(i+1,ind_lb) = sin((pi/2)*meyeraux((freq(ind_lb)-lb)/(alpha^(i-1)*df)));
        ind_ub = find((freq>=f_middle(i))&(freq<=ub));
        filter_bank(i+1,ind_ub) = cos((pi/2)*meyeraux((freq(ind_ub)-f_middle(i))/(alpha^i*df)));
    end

    % last band
    lb = f_middle(end); ub = fh+transit_width;
    ind_lb = find((freq>=lb)&(freq<=fh));
    filter_bank(end,ind_lb) = sin((pi/2)*meyeraux((freq(ind_lb)-lb)/(fh-lb)));
    ind_ub = find((freq>=fh)&(freq<=ub));
    filter_bank(end,ind_ub) = cos((pi/2)*meyeraux((freq(ind_ub)-fh)/transit_width));

    % compute fft
    gdat_fft = fft(gdat, L_fft);
    gdat_rfft = gdat_fft(1:L_fft/2+1);

    % filter by multiply by each window
    gdat_filtered = repmat(gdat_rfft,n_bands+2,1) .* filter_bank;
    gdat_filtered(:,1) = gdat_filtered(:,1)./2.0;
    gdat_filtered(:,end) = gdat_filtered(:,end)./2.0;

    gdat_filtered_fullfft = [2.0*gdat_filtered, zeros(n_bands+2, L_fft-size(gdat_filtered,2))];

    gdat_hilbert = ifft(gdat_filtered_fullfft, [], 2);
    gdat_hilbert_orig = gdat_hilbert(:,1:size(gdat,2));

    gdat_band = abs(gdat_hilbert_orig).^2;
    gdat_band = gdat_band ./ mean(gdat_band,2);
    filtered_signal = sqrt(mean(gdat_band, 1));
end
