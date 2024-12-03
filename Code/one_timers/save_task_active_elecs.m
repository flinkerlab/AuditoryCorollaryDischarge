clc; clear; close all;

% add path to visualization-tools
vistools_path = '/Volumes/research/Epilepsy_ECOG/Personal/KhalilianAmir/developer/visualization-tools/matlab';
addpath(genpath(vistools_path))

% set up the subject
root_path = pwd;
Subj = 'NY708';
tasks2consider = 'AudRep';
preprocessing_mode = 'MutliBandHG';
number_log_dist_bands = 8;
flag_clean_trials = false;
flag_auto_remove_from_time_zscore = true;

colors = get_per_task_colors(); 

fprintf('working on tasks:\n');
fprintf([tasks2consider, '\n']);

EP = EcogPreProcessing(root_path, Subj, tasks2consider,...
                       'epoch_start',-512,...
                       'epoch_end', 513,...
                       'flag_response_lock', true,...
                       'preprocessing_mode', preprocessing_mode,...
                       'number_log_dist_bands', number_log_dist_bands,...
                       'flag_clean_trials', flag_clean_trials,...
                       'Areas2Consider', {}, ...
                       'flag_auto_remove_from_time_zscore', flag_auto_remove_from_time_zscore);
[EP, ecogdata, stamps, areas, coords] = EP.get_ecogdata();

% set the figure
SCR_SZ = get(0, 'Screensize');
fig = figure('Position',SCR_SZ);
nrows = ceil(sqrt(size(ecogdata,1)));
ncols = nrows;
tiledlayout(nrows,ncols, 'Padding', 'none', 'TileSpacing', 'compact');

% loop over elecs
t = linspace(EP.epoch_start, EP.epoch_end, stamps.time_length)/EP.srate;
bad_elecs = stamps.glob.bad_elecs;
active_elecs = EcogPreProcessing.get_active_elec(ecogdata,...
                                  1.0, stamps.glob.bad_elecs);
fprintf('number of active elec: %2d\n', length(active_elecs));
for k = 1:size(ecogdata,1)
    nexttile(k);
    stdshade(squeeze(ecogdata(k,:,:)),...
             0.35, colors.(tasks2consider), t, 0);
    ylim([-0.5,1.5]);
    hold on;
    if ismember(k, bad_elecs); 
        elec_color='red';
    elseif ismember(k, active_elecs); 
        elec_color= colors.(tasks2consider);
    else; 
        elec_color='black'; 
    end
    title(sprintf('%3d,%s, %s',...
                  k, stamps.elec_labels{k}, areas{k}), 'Color', elec_color);
end

% plot selected elec on the brain
% visualization files
BrainFile = './Data/VisualData/fsaverage/lh.pial';
AnnotFile = './Data/VisualData/fsaverage/lh.aparc.split_STG_MTG.annot';
VT = visualtools('Subj', 'fs-average', 'HS', 'lh',...
                 'BrainFile', BrainFile,...
                 'AnnotFile', AnnotFile);
% convert coords to fs
coords.fs_coords = T1_to_fsaverage_utils(Subj, coords);

elec_color = ones(size(ecogdata,1),3);
elec_color(active_elecs,:) = repmat(colors.(tasks2consider),length(active_elecs),1);

VT.PlotElecOnBrain(coords.fs_coords, 'ElecColor', elec_color);

% save active electrodes 
% prepare active electrode inds to be saved
fn_save = fullfile(stamps.glob.SJdir, 'analysis', 'Active_elecs.mat');
SE.active_elecs = stamps.selected_elecs_orig_ind(active_elecs);
SE.active_elecs_labels = stamps.elec_labels(active_elecs);
if exist(fn_save, 'file')
    warning(sprintf('file %s already exists\n making backup copy\n',fn_save));
    fn_backup = sprintf('%s_backup.mat', fn_save(1:end-4));
    copyfile(fn_save, fn_backup);
end
save(fn_save, '-struct', 'SE');

% test the saved file wrt Areas2Consider

EP_ = EcogPreProcessing(root_path, Subj, tasks2consider,...
                        'epoch_start',-512,...
                        'epoch_end', 513,...
                        'flag_response_lock', true,...
                        'preprocessing_mode', preprocessing_mode,...
                        'number_log_dist_bands', number_log_dist_bands,...
                        'flag_clean_trials', flag_clean_trials,...
                        'flag_auto_remove_from_time_zscore', flag_auto_remove_from_time_zscore,...
                        'flag_load_only_active_elecs', true);
[EP_, ecogdata_, stamps_, areas_, coords_] = EP_.get_ecogdata();

% plot signals  
% set the figure
SCR_SZ = get(0, 'Screensize');
fig = figure('Position',SCR_SZ);
nrows = ceil(sqrt(size(ecogdata_,1)));
ncols = nrows;
tiledlayout(nrows,ncols, 'Padding', 'none', 'TileSpacing', 'compact');

% loop over elecs
t = linspace(EP_.epoch_start, EP_.epoch_end, stamps_.time_length)/EP_.srate;
for k = 1:size(ecogdata_,1)
    nexttile(k);
    stdshade(squeeze(ecogdata_(k,:,:)),...
             0.35, colors.(tasks2consider), t, 0);
    ylim([-0.5,1.5]);
    title(sprintf('%3d,%s, %s',...
                  stamps_.selected_elecs_orig_ind(k), stamps_.elec_labels{k}, areas_{k}), 'Color', 'black');
end

% plot selected elec on the brain
% visualization files
BrainFile = './Data/VisualData/fsaverage/lh.pial';
AnnotFile = './Data/VisualData/fsaverage/lh.aparc.split_STG_MTG.annot';
VT = visualtools('Subj', 'fs-average', 'HS', 'lh',...
                 'BrainFile', BrainFile,...
                 'AnnotFile', AnnotFile);
% convert coords to fs
coords_.fs_coords = T1_to_fsaverage_utils(Subj, coords_);

elec_color = repmat(colors.(tasks2consider),size(ecogdata_,1),1);

VT.PlotElecOnBrain(coords_.fs_coords, 'ElecColor', elec_color);
