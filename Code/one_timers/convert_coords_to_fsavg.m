clc; clear; close all;

% add path to visualization-tools
vistools_path = '/Volumes/research/Epilepsy_ECOG/Personal/KhalilianAmir/developer/visualization-tools/matlab';
addpath(genpath(vistools_path))

% set up the subject
root_path = pwd;
Subj = 'NY749';

csv_dir = fullfile(root_path, 'Data', Subj,...
                   'analysis', 'coordinates.csv');

coordinates = readtable(csv_dir);
coords.MNI = table2array(coordinates(:,2:4));
coords.T1 = table2array(coordinates(:,6:8));
coords.labels = coordinates{:,1};

fs_coords = T1_to_fsaverage(Subj, coords);

coordinates.fs_coords_1 = fs_coords(:,1);
coordinates.fs_coords_2 = fs_coords(:,2);
coordinates.fs_coords_3 = fs_coords(:,3);


% csv_wrt = fullfile(root_path, 'Data', Subj,...
%                    'analysis', 'coordinates_fs.csv');

% writetable(coordinates, csv_wrt);


% set up the subject
root_path = pwd;
Subj2 = 'NY742';

path_subj = '/Volumes/research/Epilepsy_ECOG/bigpurple_autorecon/fs/autoNY742/surf';

csv_dir2 = fullfile(root_path, 'Data', Subj2,...
                   'analysis', 'coordinates.csv');

coordinates2 = readtable(csv_dir2);
coords2.MNI = table2array(coordinates2(:,2:4));
coords2.T1 = table2array(coordinates2(:,6:8));
coords2.labels = coordinates2{:,1};

fs_coords2 = T1_to_fsaverage(Subj2, coords2, 'path_subj', path_subj);

coordinates2.fs_coords_1 = fs_coords2(:,1);
coordinates2.fs_coords_2 = fs_coords2(:,2);
coordinates2.fs_coords_3 = fs_coords2(:,3);


