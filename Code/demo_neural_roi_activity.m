clc; clear; close all;

% set up the subject
root_path = pwd;
Subj = 'NY742';
tasks2consider = 'AudRep';
preprocessing_mode = 'ChangefromBaseline';
number_log_dist_bands = 1;
flag_response_lock = true;
flag_clean_trials = false;
selected_elecs = [7,15,23,22,30,38,102,103,111,113,118,...
                  120,125,54,45,116,117,44,51,124,128,...
                  19,27,28,37,74,95,98,108,109,112,20,29,...
                  65,66,70,71,72,75,99,101,110,18,92,67,104];
RoIs2Plot = {'STG', 'IFG', 'precentral', 'postcentral'};

fprintf('working on tasks:\n');
fprintf([tasks2consider, '\n']);

EP = EcogPreProcessing(root_path, Subj, tasks2consider,...
                       'epoch_start',-512,...
                       'epoch_end', 512,...
                       'flag_response_lock', flag_response_lock,...
                       'preprocessing_mode', preprocessing_mode,...
                       'number_log_dist_bands', number_log_dist_bands,...
                       'flag_clean_trials', flag_clean_trials,...
                       'Areas2Consider', {}, ...
                       'selected_elecs', selected_elecs);
[EP, ecogdata, stamps, coords] = EP.get_ecogdata();

% set the figure
SCR_SZ = get(0, 'Screensize');
SCR_SZ(end-1) = SCR_SZ(end-1)/2;
SCR_SZ(end) = SCR_SZ(end)/3;
fig = figure('Position',SCR_SZ);
tiledlayout(1, 1, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile();

RoIsColor = turbo(length(RoIs2Plot));

% get elecs in each area
inds_in_areas = EP.get_elecs_in_area(coords.areas);

for i = 1:length(RoIs2Plot)
    roi_signal = ecogdata(inds_in_areas.(RoIs2Plot{i}),:,:);
    roi_signal = squeeze(mean(roi_signal,2));
    hold on;
    [line,fill] = stdshade(roi_signal,...
                        0.35, RoIsColor(i,:),...
                        stamps.time_vector, 20);
    line.DisplayName = RoIs2Plot{i};
    fill.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
grid on;
legend('FontSize',14);
set(gca, 'FontSize', 14);
xlabel('Time from articulation (sec)');
ylabel('Ratio change from baseline')
