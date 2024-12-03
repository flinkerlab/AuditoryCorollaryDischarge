function [X, stamps, glob] = load_preproc_elecs(obj, glob)

% get gdat from gdat_CAR.mat
if strcmp(obj.reference_mode, 'gdat_CAR')
    if exist([glob.DTdir, filesep, 'gdat_CAR.mat'], 'file')
        load([glob.DTdir, filesep, 'gdat_CAR.mat']);
    else
        disp(['file requested: ', glob.DTdir, filesep, 'gdat_CAR.mat']);
        error('file not found!');
    end
elseif strcmp(obj.reference_mode, 'gdat')
    if exist([glob.DTdir, filesep, 'gdat.mat'], 'file')
        load([glob.DTdir, filesep, 'gdat.mat']);
    else
        disp(['file requested: ', glob.DTdir, filesep, 'gdat.mat']);
        error('file not found!');
    end
else
    error('reference_mode can be gdat or gdat_CAR')
end

% load Events file
if exist(glob.Eventsdir, 'file')
	load(glob.Eventsdir);
else
	disp(['file requested: ', glob.Eventsdir]);
	error('file not found!');
end

% get onset/offset 
if isfield(glob, 'ANsrate')
    fprintf('loading orig ANsrate events\n')
    if ismember(obj.task, {'AudRep', 'PicN', 'AudN', 'SenComp', 'VisRead', 'AudRep_passive'})
        % we might need onset or onset_r
        onset = double(round(obj.Qevents(Events,glob.task,{'onset'})));
        offset = double(round(obj.Qevents(Events,glob.task,{'offset'})));
        if isfield(Events, 'onset_r')
            onset_r = double(round(obj.Qevents(Events,glob.task,{'onset_r'})));
        else
            onset_r = nan(size(onset));
            if obj.flag_response_lock
                error('Asked for locked to response but no onset_r detected!');
            end
        end
        if isfield(Events, 'offset_r')
            offset_r = double(round(obj.Qevents(Events,glob.task,{'offset_r'})));
        else
            offset_r = nan(size(onset));
        end
        if obj.flag_response_lock
            % we are locked to speech onset
            stamps.resp_on = onset_r + obj.epoch_start;
            if strcmp(obj.epoch_end, 'auto');
                trial_time = offset_r - onset_r;
                if ~all(trial_time>=0);
                    error('offset_r - onset_r < 0 found!')
                end
                time_length = ceil(mean(trial_time));
                offset_r = onset_r + time_length;
                stamps.resp_off = offset_r;
            elseif isnumeric(obj.epoch_end)
                time_length = obj.epoch_end - obj.epoch_start;
                offset_r = onset_r + obj.epoch_end;
                stamps.resp_off = offset_r;
            else
                error('I do not know what you want');
            end
            stamps.stim_on = onset;
            stamps.stim_off = offset;
        else
            % we are locked to stim 
            stamps.stim_on = onset + obj.epoch_start;
            if strcmp(obj.epoch_end, 'auto');
                trial_time = offset - onset;
                if ~all(trial_time>=0);
                    error('offset - onset < 0 found!')
                end
                time_length = ceil(mean(trial_time));
                offset = onset + time_length;
                stamps.stim_off = offset;
            elseif isnumeric(obj.epoch_end)
                time_length = obj.epoch_end - obj.epoch_start;
                offset = onset + obj.epoch_end;
                stamps.stim_off = offset;
            else
                error('I do not know what you want');
            end
            stamps.resp_on = onset_r;
            stamps.resp_off = offset_r;
        end
        % set baseline times
        stamps.baseline_on = onset + obj.baseline_start;
        stamps.baseline_off = onset + obj.baseline_end;
    else
        error(sprintf('task: %s is not in valid list', obj.task));
    end
else
	error('srate not found in ANsrate not found in glob!');
end

% unpack params
% filter frequncy bounds
f1 = obj.frequency_band_low;
f2 = obj.frequency_band_high;
num_freq = obj.number_log_dist_bands;


% initialize matrix X
badevents = obj.Qevents(Events,glob.task,{'badevent'}); 
stamps.time_length = time_length+1;
stamps.time_vector = linspace(obj.epoch_start, obj.epoch_end, stamps.time_length)/obj.srate;
X = zeros(size(gdat,1), length(stamps.resp_on)-sum(badevents), time_length+1);

% loop over requested electrods
if obj.flag_verbose; 
	fprintf('pre-processing the data... mode:%s\n', obj.preprocessing_mode);
	textprogressbar('Electrodes processed: ');
end

% find times that are out of gamma dis
if obj.flag_auto_remove_from_time_zscore
    [time_inds_remove, thresh] = obj.find_outof_gamma_dist(gdat, glob.srate);
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
            signal_in_mid = abs(obj.hilbert_filter(band,glob.srate,...
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
        evti = 1;
        for i = 1:length(stamps.resp_on)
            if badevents(i); continue; end
            if obj.flag_response_lock
                X(j, evti, :) = band_avg(stamps.resp_on(i):stamps.resp_off(i));
                evti = evti+1;
            else
                X(j, evti, :) = band_avg(stamps.stim_on(i):stamps.stim_off(i));
                evti = evti+1;
            end
        end
    elseif strcmp(obj.preprocessing_mode, 'MutliBandHGfromBaseline');
        freq_bands = logspace(log10(f1),log10(f2),num_freq+1);
        band_avg = zeros(size(band));
        for frq_i = 1:length(freq_bands)-1
            band_filtered = abs(obj.hilbert_filter(band,glob.srate,...
                                      freq_bands(frq_i),freq_bands(frq_i+1)));
            band_avg = band_avg + band_filtered./mean(band_filtered);
        end
        band_avg = band_avg ./ num_freq;
        evti = 1;
        for i = 1:length(stamps.resp_on)
            if badevents(i); continue; end
            baseline_band = band_avg(stamps.baseline_on(i):stamps.baseline_off(i));
            if obj.flag_response_lock
                X(j, evti, :) = (band_avg(stamps.resp_on(i):stamps.resp_off(i)) - mean(baseline_band,2))./mean(baseline_band,2);
                evti = evti+1;
            else
                X(j, evti, :) = (band_avg(stamps.stim_on(i):stamps.stim_off(i)) - mean(baseline_band,2))./mean(baseline_band,2);
                evti = evti+1;
            end
        end
    elseif strcmp(obj.preprocessing_mode, 'MeyerMultiBand');
        band_avg = obj.meyer_multiband(band, glob.srate, f1, f2, num_freq);
        evti = 1;
        for i = 1:length(stamps.resp_on)
            if badevents(i); continue; end
            baseline_band = band_avg(stamps.baseline_on(i):stamps.baseline_off(i));
            if obj.flag_response_lock
                X(j, evti, :) = (band_avg(stamps.resp_on(i):stamps.resp_off(i)) - mean(baseline_band,2))./mean(baseline_band,2);
                evti = evti+1;
            else
                X(j, evti, :) = (band_avg(stamps.stim_on(i):stamps.stim_off(i)) - mean(baseline_band,2))./mean(baseline_band,2);
                evti = evti+1;
            end
        end
    else
        band =  abs(obj.hilbert_filter(band,glob.srate,f1,f2));
    end
    if strcmp(obj.preprocessing_mode, 'ZScoreAll');
        band = zscore(band);
        evti = 1;
        for i = 1:length(stamps.resp_on)
            if badevents(i); continue; end
            if obj.flag_response_lock
                X(j, evti, :) = band(stamps.resp_on(i):stamps.resp_off(i));
                evti = evti+1;
            else
                X(j, evti, :) = band(stamps.stim_on(i):stamps.stim_off(i));
                evti = evti+1;
            end
        end
    elseif strcmp(obj.preprocessing_mode, 'ZScoreTrials');
        evti = 1;
        for i = 1:length(stamps.resp_on)
            if badevents(i); continue; end
            if obj.flag_response_lock
                X(j, evti, :) = band(stamps.resp_on(i):stamps.resp_off(i));
                evti = evti+1;
            else
                X(j, evti, :) = band(stamps.stim_on(i):stamps.stim_off(i));
                evti = evti+1;
            end
        end
        % z-score after trial structure is created;
        X(j,:,:) = zscore(X(j,:,:), 0, 'all');
    elseif strcmp(obj.preprocessing_mode, 'GammaNormal');
        % fit a gamma dist to the entire band and devide by mean
        [p,c] = gamfit(band);
        band = gamcdf(band,p(1),p(2));
        band = norminv(band);
        evti = 1;
        for i = 1:length(stamps.resp_on)
            if badevents(i); continue; end
            if obj.flag_response_lock
                X(j, evti, :) = band(stamps.resp_on(i):stamps.resp_off(i));
                evti = evti+1;
            else
                X(j, evti, :) = band(stamps.stim_on(i):stamps.stim_off(i));
                evti = evti+1;
            end
        end
    elseif strcmp(obj.preprocessing_mode, 'ChangefromBaseline')
        evti = 1;
        for i = 1:length(stamps.resp_on)
            if badevents(i); continue; end
            baseline_band = band(stamps.baseline_on(i):stamps.baseline_off(i));
            baseline_band_mean = mean(baseline_band,2);
            if obj.flag_response_lock
                X(j, evti, :) = (band(stamps.resp_on(i):stamps.resp_off(i)) - baseline_band_mean)./baseline_band_mean;
                evti = evti+1;
            else
                X(j, evti, :) = (band(stamps.stim_on(i):stamps.stim_off(i)) - baseline_band_mean)./baseline_band_mean;
                evti = evti+1;
            end
        end
    elseif strcmp(obj.preprocessing_mode, 'Regular');
        evti = 1;
        for i = 1:length(stamps.resp_on)
            if badevents(i); continue; end
            if obj.flag_response_lock
                X(j, evti, :) = band(stamps.resp_on(i):stamps.resp_off(i));
                evti = evti+1;
            else
                X(j, evti, :) = band(stamps.stim_on(i):stamps.stim_off(i));
                evti = evti+1;
            end
        end
    elseif strcmp(obj.preprocessing_mode, 'MutliBandHG') || strcmp(obj.preprocessing_mode, 'MutliBandHGfromBaseline')|| strcmp(obj.preprocessing_mode, 'MeyerMultiBand')
    else
        error(sprintf('pre-processing mode %s not valid', obj.preprocessing_mode))
    end
    if obj.flag_verbose
        textprogressbar(j/size(gdat,1)*100);
    end
end
        
if obj.flag_verbose
	textprogressbar('done.');
end

end
