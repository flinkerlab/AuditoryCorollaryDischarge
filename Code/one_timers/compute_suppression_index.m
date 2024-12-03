function [SI_mean, SI_se] = compute_suppression_index(Subj, selected_elecs, varargin)
    p = inputParser;
    addParameter(p, 'root_path', pwd);
    addParameter(p, 'tasks2consider', 'AudRep');
    addParameter(p, 'preprocessing_mode', 'Regular');
    parse(p, varargin{:});
    root_path = p.Results.root_path;
    tasks2consider = p.Results.tasks2consider;
    preprocessing_mode = p.Results.preprocessing_mode;

    % load the reference signals (i.e. hearing)
    EP_ref = EcogPreProcessing(root_path, Subj, tasks2consider,...
                               'epoch_start', 0,...
                               'epoch_end', round(0.3*512),...
                               'flag_response_lock', false,...
                               'preprocessing_mode', preprocessing_mode,...
                               'Areas2Consider', {}, ...
                               'selected_elecs', selected_elecs);
    [EP_ref, ecogdata_ref, stamps_ref, coords_ref] = EP_ref.get_ecogdata();

    % load the task signals (i.e. speaking)
    EP_tsk = EcogPreProcessing(root_path, Subj, tasks2consider,...
                               'epoch_start', 0,...
                               'epoch_end', round(0.3*512),...
                               'flag_response_lock', true,...
                               'preprocessing_mode', preprocessing_mode,...
                               'Areas2Consider', {}, ...
                               'selected_elecs', selected_elecs);
    [EP_tsk, ecogdata_tsk, stamps_tsk, coords_tsk] = EP_tsk.get_ecogdata();

    % compute the mean activity per electrode per trial
    ref_mean = squeeze(mean(ecogdata_ref, 3));
    tsk_mean = squeeze(mean(ecogdata_tsk, 3));

    % compute the suppression index (SI) per trial per electrode
    SI = (ref_mean - tsk_mean)./(ref_mean + tsk_mean);
    % compute the mean and SEM for SI over trials
    SI_mean = squeeze(mean(SI, 2));
    SI_se = squeeze(std(SI,[],2))./sqrt(size(SI,2));
end
