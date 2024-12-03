clc; clear; close all;

subj = 'NY717';
task = 'PicN'; 

fn_crct = fullfile('./Data', subj, 'analysis', task, 'Events.mat');
Ecrct = load(fn_crct, 'Events');
Ecrct = Ecrct.Events;

fn_old = fullfile('/Users/khalia03/Documents/ECoGAR_revision/Data', subj, 'analysis', task, 'Events.mat');
Eold = load(fn_old, 'Events');
Eold = Eold.Events;

onset_r_crct = [Ecrct(:).onset_r];
onset_r_old = double([Eold(:).onset_r]);

plot(onset_r_crct-onset_r_old, 'k-*', 'MarkerSize', 8);
