clc; clear; close all;

% set up the subject
root_path = pwd;
Subj = 'NY742';
tasks2consider = 'AudRep';
preprocessing_mode = 'Regular';
number_log_dist_bands = 1;
selected_elecs = [7,15,23,22,30,38,102,103,111,113,118,120,125,54,45,116,117,44,51,124,128,19,27,28,37,74,95,98,108,109,112,20,29,65,66,70,71,72,75,99,101,110,18,92,67,104];
% select only STG elecs
selected_elecs = selected_elecs(1:14);

fprintf('working on tasks:\n');
fprintf([tasks2consider, '\n']);

% load the reference signals (i.e. hearing)
EP_ref = EcogPreProcessing(root_path, Subj, tasks2consider,...
                           'epoch_start', -round(9/200*512),...
                           'epoch_end', round(0.3*512)+round(10/200*512),...
                           'flag_response_lock', false,...
                           'preprocessing_mode', preprocessing_mode,...
                           'number_log_dist_bands', number_log_dist_bands,...
                           'Areas2Consider', {}, ...
                           'flag_resample_data', true,...
                           'resample_srate', 200,...
                           'selected_elecs', selected_elecs);
[EP_ref, ecogdata_ref, stamps_ref, coords_ref] = EP_ref.get_ecogdata();

% load the task signals (i.e. speaking)
EP_tsk = EcogPreProcessing(root_path, Subj, tasks2consider,...
                           'epoch_start', -round(9/200*512),...
                           'epoch_end', round(0.3*512)+round(10/200*512),...
                           'flag_response_lock', true,...
                           'preprocessing_mode', preprocessing_mode,...
                           'number_log_dist_bands', number_log_dist_bands,...
                           'Areas2Consider', {}, ...
                           'flag_resample_data', true,...
                           'resample_srate', 200,...
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
figure;
[SI_sorted, sort_ind] = sort(SI_mean);
plot(1:size(SI_sorted), SI_sorted, 'k*', 'MarkerSize', 9);
hold on;
errorbar(1:size(SI_sorted), SI_sorted, SI_se(sort_ind),...
         'k.', 'LineWidth', 2);

% set the figure
SCR_SZ = get(0, 'Screensize');
fig = figure('Position',SCR_SZ);
nrows = ceil(sqrt(size(ecogdata_ref,1)));
ncols = nrows;
tiledlayout(nrows,ncols, 'Padding', 'none', 'TileSpacing', 'compact');

% loop over elecs and plot the mean activity for each electrode
bad_elecs = stamps_ref.glob.bad_elecs;
for k = 1:size(ecogdata_ref,1)
    nexttile(k);
    stdshade(squeeze(ecogdata_ref(k,:,:)),...
             0.35, [0.5 0.5 0.5], stamps_ref.time_vector, 0);
    hold on;
    stdshade(squeeze(ecogdata_tsk(k,:,:)),...
             0.35, [0.5 0.5 1.0], stamps_tsk.time_vector, 0);
    if ismember(EP_ref.selected_elecs(k), bad_elecs); 
        elec_color='red';
    else; 
        elec_color='black'; 
    end
    title(sprintf('%3d, %s, %s, SI: %2.2f',...
                  k, coords_ref.labels{k},...
                  coords_ref.areas{k},SI_mean(k)),...
                  'Color', elec_color);
end


