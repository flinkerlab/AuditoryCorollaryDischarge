classdef EcogPreProcessing
	% properties for the EcogPreProcessing class
	properties
		% properties for data
		root_path
		subj
		task
		% properties for pre-processing
        selected_elecs = [];
		epoch_start = -10;
		epoch_end = 512;
        baseline_start = -160;
        baseline_end = -32;
		srate = 512;
		flag_response_lock = true;
		frequency_band_low = 70;
		frequency_band_high = 150;
        number_log_dist_bands = 8;
		preprocessing_mode = 'ZScoreAll'; % z-score all or z-score trial 'ZScoreTrails' 'Regular'
        reference_mode = 'gdat_CAR';
        flag_auto_remove_from_time_zscore = false;
        flag_verbose = true;
        flag_clean_trials = false;
        flag_use_run_2 = false;
        flag_load_only_active_elecs = false;
        flag_resample_data = false;
        resample_srate = 512;
        clean_trials_thresh = 15;
        Areas2Consider = {'cSTG','rSTG','mSTG',...
                          'cMTG','rMTG','mMTG',...
                          'caudalmiddlefrontal',...
                          'insula',...
                          'parsopercularis',...
                          'parsorbitalis',...
                          'parstriangularis',...
                          'precentral',...
                          'postcentral',...
                          'rostralmiddlefrontal',...
                          'supramarginal'};
		% ecog data
		ecogdata = [];
		stamps = struct();
	end

	% methods for the EcogPreProcessing class
	methods
		% Constructor class
		function [obj] = EcogPreProcessing(root_path, subj, task, varargin)
			% method to construct the EcogPreProcessing class
			% parse the input to the struct
			p = inputParser;
			p.KeepUnmatched = true;
			ValidTasks = {'AudRep', 'PicN', 'AudN', 'SenComp', 'VisRead', 'AudRep_passive', 'AudRep_imagine'}; 
			addRequired(p, 'root_path', @(s) (isdir(s)));
			addRequired(p, 'subj', @(s) (ischar(s))||(isstring(s)&&length(s)==1));
			addRequired(p, 'task', @(s) any(validatestring(s, ValidTasks)));
            addParameter(p, 'selected_elecs', obj.selected_elecs);
			addParameter(p, 'epoch_start',  obj.epoch_start);
			addParameter(p, 'epoch_end', obj.epoch_end);
            addParameter(p, 'baseline_start', obj.baseline_start);
            addParameter(p, 'baseline_end', obj.baseline_end);
			addParameter(p, 'srate', obj.srate);
			addParameter(p, 'flag_response_lock', obj.flag_response_lock);
			addParameter(p, 'frequency_band_low', obj.frequency_band_low);
			addParameter(p, 'frequency_band_high', obj.frequency_band_high);
			addParameter(p, 'number_log_dist_bands', obj.number_log_dist_bands);
			addParameter(p, 'preprocessing_mode', obj.preprocessing_mode);
			addParameter(p, 'reference_mode', obj.reference_mode);
			addParameter(p, 'flag_verbose', obj.flag_verbose);
			addParameter(p, 'flag_clean_trials', obj.flag_clean_trials);
			addParameter(p, 'clean_trials_thresh', obj.clean_trials_thresh);
			addParameter(p, 'flag_auto_remove_from_time_zscore', obj.flag_auto_remove_from_time_zscore);
			addParameter(p, 'flag_load_only_active_elecs', obj.flag_load_only_active_elecs);
			addParameter(p, 'flag_resample_data', obj.flag_resample_data);
			addParameter(p, 'resample_srate', obj.resample_srate);
			addParameter(p, 'Areas2Consider', obj.Areas2Consider);
			% parse
			parse(p, root_path, subj, task, varargin{:});
			% set parsed to obj
			obj.root_path = p.Results.root_path;
			obj.subj = p.Results.subj;
			obj.task = p.Results.task;
			obj.selected_elecs = p.Results.selected_elecs;
			obj.epoch_start = p.Results.epoch_start;
			obj.epoch_end = p.Results.epoch_end;
			obj.baseline_start = p.Results.baseline_start;
			obj.baseline_end = p.Results.baseline_end;
			obj.srate = p.Results.srate;
			obj.flag_response_lock = p.Results.flag_response_lock;
			obj.frequency_band_low = p.Results.frequency_band_low;
			obj.frequency_band_high = p.Results.frequency_band_high;
			obj.number_log_dist_bands = p.Results.number_log_dist_bands;
			obj.preprocessing_mode = p.Results.preprocessing_mode;
			obj.reference_mode = p.Results.reference_mode;
			obj.flag_verbose = p.Results.flag_verbose;
			obj.flag_clean_trials = p.Results.flag_clean_trials;
            obj.clean_trials_thresh = p.Results.clean_trials_thresh;
            obj.flag_auto_remove_from_time_zscore = p.Results.flag_auto_remove_from_time_zscore;
            obj.flag_load_only_active_elecs = p.Results.flag_load_only_active_elecs;
            obj.flag_resample_data = p.Results.flag_resample_data;
            obj.resample_srate = p.Results.resample_srate;
            obj.Areas2Consider = p.Results.Areas2Consider;
		end

		% load and pre-process ecog data functions
		function [obj, ecogdata, stamps, coords] = get_ecogdata(obj)
			% this function loads the ecog data for the subject and task
			% in obj.subj, obj.task
			glob = obj.get_subj_global();
			[ecogdata, stamps, glob] = obj.load_preproc_elecs(glob);
			coordinates = readtable(glob.CSVdir);
			areas = coordinates{:,11};
            if ~isempty(obj.Areas2Consider)
                SE_in_RoI = EcogPreProcessing.get_elecs_in_areas2consider(areas,...
                                              obj.Areas2Consider, glob.bad_elecs);
            else
                SE_in_RoI = 1:size(ecogdata,1);
            end
            if obj.flag_load_only_active_elecs
                if isempty(glob.selected_elecs)
                    error('no active elec file loaded in globals');
                else
                    SE = intersect(glob.selected_elecs.active_elecs,...
                                   SE_in_RoI);
                end
            else
                SE = SE_in_RoI;
            end
            if ~isempty(obj.selected_elecs)
                SE = obj.selected_elecs;
            end
            fprintf('number of elecs selected: %d\n', length(SE));
            % select only electrodes from list SE
            ecogdata = ecogdata(SE,:,:);
            % set electrode related information to coords
            coords.MNI = table2array(coordinates(SE,2:4));
            coords.T1 = table2array(coordinates(SE,6:8));
			coords.labels = coordinates{SE,1};
            coords.areas = coordinates{SE,11};
            coords.selected_elecs_orig_ind = SE;
            % form add glob to stamps for book-keeping
			stamps.glob = glob;
            if obj.flag_clean_trials
                [good_trial, bad_trial] = obj.clean_trials(ecogdata,...
                                            obj.clean_trials_thresh);
                ecogdata(:,bad_trial,:) = [];
                stamps.good_trial_ind = good_trial;
                stamps.bad_trials_ind_removed = bad_trial;
            end
            % resample the data if requested
            if obj.flag_resample_data
                [obj, ecogdata, stamps] = obj.resample_data(ecogdata, stamps);
            end
		end

        % load the gdat_CAR only without preprocessing and epoch structure
        function [gdat] = get_gdat(obj)
            % load subject glob
			glob = obj.get_subj_global();
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
        end

        % simple pre-processing routine for a band without trial event
        [gdat] = band_preproc(obj, gdat, srate);

		% subject globals
		[glob] = get_subj_global(obj)

		% load the data and pre-process
		[ecogdata, stamps, glob] = load_preproc_elecs(obj, glob);

        % load and structure the zoom file
        [zoom_tr] = load_preproc_zoom(obj, glob, stamps, varargin)

        % resample the ecogdata into new srate
        function [obj, ecogdata_rs, stamps_rs] = resample_data(obj, ecogdata, stamps)
            % function to convert the ecogdata from fs to fs_new
            % ecogdata in shape elec x trial x time
            fs = obj.srate;
            fs_new = obj.resample_srate;
            [p,q] = rat(fs_new/fs);
            stamps_rs = struct();
            [ecogdata_rs, stamps_rs.time_vector] = resample(ecogdata,...
                                                            stamps.time_vector,...
                                                            fs_new, p, q,...
                                                            'Dimension', 3);
            % covert stamps from fs to fs_new
            stamps_rs.stim_on  = round((stamps.stim_on/q)*p);
            stamps_rs.stim_off = round((stamps.stim_off/q)*p);
            stamps_rs.resp_on  = round((stamps.resp_on/q)*p);
            stamps_rs.resp_off = round((stamps.resp_off/q)*p);
            stamps_rs.baseline_on  = round((stamps.baseline_on/q)*p);
            stamps_rs.baseline_off = round((stamps.baseline_off/q)*p);
            stamps_rs.time_length = length(stamps_rs.time_vector);
            stamps_rs.glob = stamps.glob;
            stamps_rs.glob.srate = fs_new;
            stamps_rs.glob.ANsrate = fs_new;
            % convert obj events
            obj.epoch_start = round((obj.epoch_start/q)*p);
            obj.epoch_end = round((obj.epoch_end/q)*p);
            if length(obj.epoch_start:obj.epoch_end)~=stamps_rs.time_length
                % the epoch coversion is not correct
                error('epoch conversion mismatch');
            end
            % remove samples for boundary effect
            ecogdata_rs = ecogdata_rs(:,:, 10:end-10);
            stamps_rs.time_vector = stamps_rs.time_vector(10:end-10);
            stamps_rs.time_length = length(stamps_rs.time_vector);
            obj.epoch_start = obj.epoch_start+9;
            obj.epoch_end = obj.epoch_end-10;
            if length(obj.epoch_start:obj.epoch_end)~=size(ecogdata_rs,3)
                % the epoch coversion is not correct
                error('epoch conversion mismatch in trunctation');
            end
            % change obj srate
            obj.srate = fs_new;
        end

	end

	methods (Static)
		% get events in a Database
		[out] = Qevents(Database, ev, fields,str_fld,str_val);

		% hilbert filtering of the signals
		[filt_signal, wind_out] = hilbert_filter(input, sampling_rate,...
									             lower_bound, upper_bound);

        % meyer filtering
        [filtered_signal] = meyer_multiband(gdat, Fs, fl, fh, n_bands)

		% band-padd filtering of the signal
		[filt_signal, Pl] = band_pass(input, sampling_rate,...
									  lower_bound, upper_bound);

        % find times that are out of dist
        [time_inds, thresh] = find_outof_gamma_dist(gdat, srate);


        % clean trials
        function [good_trial, bad_trial] = clean_trials(ecogdata, thresh)
            ind = find(ecogdata>=thresh);
            [i1,i2,i3] = ind2sub(size(ecogdata), ind);
            i2 = unique(i2);
            good_trial= setdiff(1:size(ecogdata,2),i2);
            bad_trial = i2;
        end

        % get electrode indices within regions of interest
        function [elec_inds] = get_elecs_in_areas2consider(elec_areas, Areas2Consider, bad_elecs)
            elec_inds = find(ismember(elec_areas, Areas2Consider));
            elec_inds = setdiff(elec_inds, bad_elecs);
        end

        % get active electrodes based on mean signal
        function [elec_inds] = get_active_elec(ecogdata, thresh, bad_elecs);
            % compute the mean over trials
            xm = squeeze(mean(ecogdata, 2));
            % smooth the mean signal with a wavelet filter
            xm_smooth = zeros(size(xm));
            for i = 1:size(xm,1)
                if ismember(i, bad_elecs);
                    continue;
                end
                xm_smooth(i,:) = WTDN(xm(i,:), 'db8', 0.5);
            end
            std_all = std(xm_smooth, [], 2);
            [rice_mean,rice_sigma] = ricefit(std_all);
            elec_inds = find(std_all>=thresh*rice_sigma);
            elec_inds = setdiff(elec_inds, bad_elecs);
        end

        function [inds_per_area] = get_elecs_in_area(areas, varargin)
            % function to get the elecs in an ROI
            p = inputParser;
            addParameter(p, 'flag_mergeTemporal', true);
            addParameter(p, 'flag_mergeIFG', true);
            parse(p, varargin{:});
            flag_mergeTemporal = p.Results.flag_mergeTemporal;
            flag_mergeIFG = p.Results.flag_mergeIFG;
            if flag_mergeTemporal
                % merge cSTG, mSTG, rSTG into STG
                areas(find(contains(areas, 'STG'))) = {'STG'};
                % merge cMTG, mMTG, rMTG into MTG
                areas(find(contains(areas, 'MTG'))) = {'MTG'};
            end
            if flag_mergeIFG
                areas(find(contains(areas, 'parsopercularis'))) = {'IFG'};
                areas(find(contains(areas, 'parstriangularis'))) = {'IFG'};
            end
            ua = unique(areas);
            inds_per_area = struct();
            for i = 1:length(ua);
                inds_per_area.(ua{i}) = find(contains(areas, ua{i}));
            end
        end
								  
	end
end
