function visual_clustering(V, Ulb, U, coords, varargin)
    % function to visualize the connectivity results
    p = inputParser;
    addParameter(p, 'clusters_to_show', []);
    addParameter(p, 'time_start_window', []);
    DEF_vis_path = fullfile(pwd,...
                            'Code/visualization_tools/matlab');
    addParameter(p, 'vistools_path', DEF_vis_path);
    DEF_BrainFile = fullfile(DEF_vis_path,...
                             'SampleData/MNI/ch2_template_lh_pial_120519.mat');
    addParameter(p, 'BrainFile', DEF_BrainFile);
    DEF_AnnotFile = fullfile(DEF_vis_path,...
                             'SampleData/MNI/ch2_template.lh.aparc.split_STG_MTG.annot');
    addParameter(p, 'AnnotFile', DEF_AnnotFile);
    addParameter(p, 'Subj', 'MNI');
    addParameter(p, 'HS', 'lh');
    addParameter(p, 'flag_UseAnnots', true);
    addParameter(p, 'flag_MergeSTG', true);
    addParameter(p, 'flag_renormalize', false);
    parse(p, varargin{:});
    clusters_to_show = p.Results.clusters_to_show;
    if isempty(clusters_to_show)
        clusters_to_show = [1:size(V,1)];
    end
    time_start_window = p.Results.time_start_window;
    if isempty(time_start_window)
        time_start_window = [1:size(V,2)];
    end
    flag_renormalize = p.Results.flag_renormalize;
    vistools_path = p.Results.vistools_path;
    BrainFile = p.Results.BrainFile;
    AnnotFile = p.Results.AnnotFile;
    Subj = p.Results.Subj;
    HS = p.Results.HS;
    flag_UseAnnots = p.Results.flag_UseAnnots;
    flag_MergeSTG = p.Results.flag_MergeSTG;
    addpath(genpath(vistools_path));
    

    % get the weights for plotting
    [Wa, ~, ~] = EcogONMF.get_clustering_weights(U,Ulb);
    if flag_renormalize
        Wa(Wa>=0) = Wa(Wa>=0) ./ max(Wa(Wa>=0),[],'all');
        Wa(Wa<=0) = Wa(Wa<=0) ./ max(abs(Wa(Wa<=0)),[],'all');
        clim = [-1,1];
    else
        clim = [];
    end

    % create a visualization object
    VT = visualtools('Subj', Subj, 'HS', HS,...
                     'flag_UseAnnots', flag_UseAnnots,...
                     'flag_MergeSTG', flag_MergeSTG,...
                     'BrainFile', BrainFile,...
                     'AnnotFile', AnnotFile);
    if strcmp(Subj, 'MNI') || strcmp(Subj,'MNI-FS')
        coords_to_use = 'MNI';
    else
        coords_to_use = 'T1';
    end
    cmap = get_b2r_cmap();

    % make a new figure;
    SCR_SZ = get(0, 'Screensize');
    SCR_SZ(end-1) = SCR_SZ(end-1)/3*2;
    fig = figure('Position',SCR_SZ);
    nrows = length(clusters_to_show);
    ncols = 2;
    tiledlayout(nrows,ncols, 'Padding', 'none', 'TileSpacing', 'compact');

    % set ylims for component plots
    V_min = min(V(clusters_to_show,:),[],'all');
    V_max = max(V(clusters_to_show,:),[],'all');

    % loop over cluster comps
    for i = 1:length(clusters_to_show)
        % cluster number
        k = clusters_to_show(i);
        % plot the temporal comp
        nexttile(2*i-1);
        plot(time_start_window,...
            V(k,:),...
            'LineWidth', 2);
        ylim([V_min, V_max]);
        if i==1;
            xlabel('Time (sec)');
            ylabel('Normalized Directed Connectivity');
        end
        % plot the weights
        nexttile(2*i);
        VT.PlotWeightedBrain(coords.(coords_to_use),...
                             Wa(:,k),...
                             'ElecRegions', coords.areas,...
                             'Sigma', 5,...
                             'Mode', 'SimpleProject',...
                             'cmap', cmap, 'clim', clim);
    end


end
