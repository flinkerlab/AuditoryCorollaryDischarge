clc; clear; close all;

% select subject
root_path = pwd;
Subj = 'NY742';
task = 'AudRep';
selected_elecs = 7;
trial = 5;

% get the data matrix 
EP = EcogPreProcessing(root_path, Subj, task);
% load the data without CAR
EP.reference_mode = 'gdat';
gdat = EP.get_gdat();
% load the data withCAR
EP.reference_mode = 'gdat_CAR';
gdat_CAR = EP.get_gdat();

% get the subject global info
glob = EP.get_subj_global();

% get the trial Events
Events = load(glob.Eventsdir);

onsets = double(round(EP.Qevents(Events.Events,glob.task,{'onset'})));
badevent = double(round(EP.Qevents(Events.Events,glob.task,{'badevent'})));
if badevent(trial); error('the trial selected is a bad trial'); end

% take an epoch of data
epoch_start = -128; % 250 msec before stim onset
epoch_end = 1536;  % 3000 msec after stim onset
tm_inds = [onsets(trial)+epoch_start:onsets(trial)+epoch_end];
time_vector = [epoch_start:epoch_end]/glob.srate*1000;

% get the data --
trial_gdat = gdat(selected_elecs,tm_inds);
trial_gdat_CAR = gdat_CAR(selected_elecs,tm_inds);

% plot the original signal
figure;
plot(time_vector, trial_gdat, 'k', 'LineWidth', 1.5);
xlim([time_vector(1),time_vector(end)]);
ylim([-150,150]);
grid on;
xlabel('Time from trial onset (msec)');
ylabel('signal (\muV)');
title('recorded signal');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% plot the signal after CAR
figure;
plot(time_vector, trial_gdat_CAR, 'k', 'LineWidth', 1.5);
xlim([time_vector(1),time_vector(end)]);
ylim([-150,150]);
grid on;
xlabel('Time from trial onset (msec)');
ylabel('signal (\muV)');
title('recorded signal after CAR');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% plot the CAR signal
figure;
plot(time_vector, trial_gdat - trial_gdat_CAR, 'k', 'LineWidth', 1.5);
xlim([time_vector(1),time_vector(end)]);
ylim([-150,150]);
grid on;
xlabel('Time from trial onset (msec)');
ylabel('CAR (\muV)');
title('common avg. ref.');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% bandpass filter the signal to 70-150 Hz
highgamma_band = EP.band_pass(gdat_CAR(selected_elecs,:),...
                              glob.srate, 70, 150);
trial_highgamma = highgamma_band(tm_inds);

% plot the filtered highgamma signal
figure;
plot(time_vector, trial_highgamma, 'k', 'LineWidth', 1.5);
xlim([time_vector(1),time_vector(end)]);
ylim([-25, 25]);
grid on;
xlabel('Time from trial onset (msec)');
ylabel('high-gamma signal (\muV)');
title('high-gamma bandpass signal');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% compute the analytic amp. of highgamma
analytic_amp = abs(EP.hilbert_filter(gdat_CAR(selected_elecs,:),...
                                       glob.srate, 70, 150));
trial_analytic_amp = analytic_amp(tm_inds);

% plot the analytic amp. highgamma signal
figure;
plot(time_vector, trial_analytic_amp, 'k', 'LineWidth', 1.5);
xlim([time_vector(1),time_vector(end)]);
ylim([-25, 25]);
grid on;
xlabel('Time from trial onset (msec)');
ylabel('analytic amp. signal (\muV)');
title('analytic amplitude of high-gamma signal');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% downsample the envelope signal
[p,q] = rat(200/glob.srate);
trial_analytic_ds = resample(trial_analytic_amp, p, q);
time_vector_ds = resample(time_vector, p, q);
                                     
% remove the boundary for edge effects
trial_analytic_ds = trial_analytic_ds(10:end-10);
time_vector_ds = time_vector_ds(10:end-10);

% plot the analytic amp. highgamma signal at 200 Hz
figure;
plot(time_vector_ds, trial_analytic_ds, 'k', 'LineWidth', 1.5);
xlim([time_vector(1),time_vector(end)]);
ylim([-25, 25]);
grid on;
xlabel('Time from trial onset (msec)');
ylabel('analytic amp. signal (\muV)');
title('analytic amplitude of high-gamma signal');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% plotting power spectrums!
wins = 512/0.1; % srate(Hz) / freq resolution (Hz/bin)
overlap = wins/2;
[P_trial_gdat, f_512] = pwelch(gdat(selected_elecs,:),...
                               wins,overlap, [0:0.1:256],512.0);
P_trial_gdat_CAR = pwelch(gdat_CAR(selected_elecs,:),...
                          wins,overlap,[0:0.1:256],512.0);
trial_CAR = gdat(selected_elecs,:) - gdat_CAR(selected_elecs,:);
P_trial_CAR = pwelch(trial_CAR,...
                     wins,overlap,[0:0.1:256],512.0);
P_trial_highgamma = pwelch(highgamma_band,...
                           wins,overlap,[0:0.1:256],512.0);
P_trial_analytic_amp = pwelch(analytic_amp,...
                              wins,overlap,[0:0.1:256],512.0);
[P_trial_analytic_ds, f_200] = pwelch(resample(analytic_amp, p, q),...
                                      wins/512*200,overlap/512*200,[0:0.1:100],200.0);

% plot the original signal spectrum
figure;
plot(f_512, pow2db(P_trial_gdat), 'k', 'LineWidth', 1.5);
xlim([0,256]);
ylim([-40,40])
grid on;
xlabel('frequency (Hz)');
ylabel('signal PSD (dB/Hz)');
title('recorded signal -- PSD');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% plot the signal after CAR
figure;
plot(f_512, pow2db(P_trial_gdat_CAR), 'k', 'LineWidth', 1.5);
xlim([0,256]);
ylim([-40,40]);
grid on;
xlabel('frequency (Hz)');
ylabel('signal PSD (dB/Hz)');
title('recorded signal after CAR -- PSD');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% plot the CAR signal
figure;
plot(f_512, pow2db(P_trial_CAR), 'k', 'LineWidth', 1.5);
xlim([0,256]);
ylim([-40,40]);
grid on;
xlabel('frequency (Hz)');
ylabel('CAR PSD (dB/Hz)');
title('common avg. ref. -- PSD');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% plot the CAR signal
figure;
plot(f_512, pow2db(P_trial_highgamma), 'k', 'LineWidth', 1.5);
xlim([0,256]);
ylim([-60,30]);
grid on;
xlabel('frequency (Hz)');
ylabel('high-gamma signal PSD (dB/Hz)');
title('high-gamma signal -- PSD');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% plot the CAR signal
figure;
plot(f_512, pow2db(P_trial_analytic_amp), 'k', 'LineWidth', 1.5);
xlim([0,256]);
ylim([-60,30]);
grid on;
xlabel('frequency (Hz)');
ylabel('high-gamma signal PSD (dB/Hz)');
title('high-gamma analytic amp -- PSD');
set(gca, 'FontName', 'Times', 'FontSize', 15);

% plot the CAR signal
figure;
plot(f_200, pow2db(P_trial_analytic_ds), 'k', 'LineWidth', 1.5);
xlim([0,256]);
ylim([-60,30]);
grid on;
xlabel('frequency (Hz)');
ylabel('high-gamma signal PSD (dB/Hz)');
title('high-gamma analytic amp -- PSD');
set(gca, 'FontName', 'Times', 'FontSize', 15);

