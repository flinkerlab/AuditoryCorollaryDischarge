function [fs_coords] = T1_to_fsaverage(subj, coords)
    % wrapper function around visualization-tools/matlab/fsaveragetools.m

    % add path to visualization-tools
    vistools_path = '/Volumes/research/Epilepsy_ECOG/Personal/KhalilianAmir/developer/visualization-tools/matlab';
    addpath(genpath(vistools_path))

    % set data path for fs-average
    data_visuals_path = '/Volumes/research/Epilepsy_ECOG/Personal/KhalilianAmir/EcogData/VisualData';
    fsaverage_path = [data_visuals_path, filesep, 'fsaverage'];

    fs_lh_pial = [fsaverage_path, filesep, 'lh.pial'];
    fs_rh_pial = [fsaverage_path, filesep, 'rh.pial'];

    fs_lh_sphere = [fsaverage_path, filesep, 'lh.sphere.reg'];
    fs_rh_sphere = [fsaverage_path, filesep, 'rh.sphere.reg'];

    % set data path for subj
    subj_path = [data_visuals_path, filesep, subj];
    subj_lh_pial = [subj_path, filesep, subj, '.lh.pial'];
    subj_rh_pial = [subj_path, filesep, subj, '.rh.pial'];

    subj_lh_sphere = [subj_path, filesep, subj, '.lh.sphere.reg'];
    subj_rh_sphere = [subj_path, filesep, subj, '.rh.sphere.reg'];

    % find index of elecs on each hemisphere
    lh_ind = find(coords.MNI(:,1)<=0);
    rh_ind = find(coords.MNI(:,1)>0);

    % find depth electrodes (assuming label starts with D)
    depth_ind = find(startsWith(coords.labels, 'd', 'IgnoreCase', true));
    surf_ind = setdiff([1:size(coords.MNI,1)], depth_ind);

    % fs-average coords place holder
    fs_coords = nan(size(coords.T1));

    % find indices on each hemisphere that are surface
    surf_lh_ind = intersect(surf_ind, lh_ind);
    surf_rh_ind = intersect(surf_ind, rh_ind);

    % create fsaveragetools object for each hemisphere
    fstools_lh = fsaveragetools('HS', 'lh',...
                                'fn_fsaverage_pial', fs_lh_pial,...
                                'fn_fsaverage_sphere', fs_lh_sphere,...
                                'fn_subj_pial', subj_lh_pial,...
                                'fn_subj_sphere', subj_lh_sphere);

    fstools_rh = fsaveragetools('HS', 'rh',...
                                'fn_fsaverage_pial', fs_rh_pial,...
                                'fn_fsaverage_sphere', fs_rh_sphere,...
                                'fn_subj_pial', subj_rh_pial,...
                                'fn_subj_sphere', subj_rh_sphere);

    % convert coords using convert_T1_to_fsaverage
    temp_coords_lh = fstools_lh.convert_T1_to_fsaverage(coords.T1(surf_lh_ind,:));
    fs_coords(surf_lh_ind,:) = temp_coords_lh;

    temp_coords_rh = fstools_rh.convert_T1_to_fsaverage(coords.T1(surf_rh_ind,:));
    fs_coords(surf_rh_ind,:) = temp_coords_rh;
end



