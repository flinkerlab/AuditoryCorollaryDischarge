function [gdat_proc] = band_preproc(obj, gdat, srate)
    % unpack params
    % filter frequncy bounds
    f1 = obj.frequency_band_low;
    f2 = obj.frequency_band_high;
    num_freq = obj.number_log_dist_bands;

    gdat_proc = nan(size(gdat));

    % loop over requested electrods
    if obj.flag_verbose; 
        fprintf('pre-processing the data... mode:%s\n', obj.preprocessing_mode);
        textprogressbar('Electrodes processed: ');
    end

    % find times that are out of gamma dis
    if obj.flag_auto_remove_from_time_zscore
        [time_inds_remove, thresh] = obj.find_outof_gamma_dist(gdat, srate);
    else
        time_inds_remove = cell(size(gdat,1),1);
        thresh = inf;
    end

    for j = 1:size(gdat,1) % looping over elecs
        % signal band for this electrode
        band = gdat(j,:); 
        if strcmp(obj.preprocessing_mode, 'MutliBandHG');
            freq_bands = logspace(log10(f1),log10(f2),num_freq+1);
            band_avg = zeros(size(band));
            for frq_i = 1:length(freq_bands)-1
                signal_in_mid = abs(obj.hilbert_filter(band,srate,...
                                        freq_bands(frq_i),freq_bands(frq_i+1)));
                if obj.flag_auto_remove_from_time_zscore
                    % remove the times in time_inds
                    temp_nan = signal_in_mid;
                    temp_nan(time_inds_remove{j}) = nan;
                    signal_in_mid = (signal_in_mid - mean(temp_nan,'omitnan'))./...
                                    std(temp_nan,[],'omitnan');
                else
                    signal_in_mid = zscore(signal_in_mid);
                end
                band_avg = band_avg + signal_in_mid;
            end
            band_avg = band_avg / num_freq;
            gdat_proc(j,:) = band_avg;
        elseif strcmp(obj.preprocessing_mode, 'ZScoreAll');
            band =  abs(obj.hilbert_filter(band,srate,f1,f2));
            band = zscore(band);
            gdat_proc(j,:) = band;
        else
            error('Invalid pre-prcosseing mode for this function')
        end
        if obj.flag_verbose
            textprogressbar(j/size(gdat,1)*100);
        end
    end
            
    if obj.flag_verbose
        textprogressbar('done.');
    end

end
