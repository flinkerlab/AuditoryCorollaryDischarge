function visual_clustering_suppression(V, Ulb, U, coords, SI, varargin)
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
    addParameter(p, 'Inflow_RoI', 'STG');
    parse(p, varargin{:});
    clusters_to_show = p.Results.clusters_to_show;
    time_start_window = p.Results.time_start_window;
    if isempty(time_start_window)
        time_start_window = [1:size(V,2)];
    end
    if isempty(clusters_to_show)
        [MaxVals, MaxInds] = max(V, [], 2);
        ind = find(time_start_window(MaxInds)<=0 & time_start_window(MaxInds)>=-0.2);
        clusters_to_show = [ind];
    end
    flag_renormalize = p.Results.flag_renormalize;
    Inflow_RoI = p.Results.Inflow_RoI;
    vistools_path = p.Results.vistools_path;
    BrainFile = p.Results.BrainFile;
    AnnotFile = p.Results.AnnotFile;
    Subj = p.Results.Subj;
    HS = p.Results.HS;
    flag_UseAnnots = p.Results.flag_UseAnnots;
    flag_MergeSTG = p.Results.flag_MergeSTG;
    addpath(genpath(vistools_path));
    

    % get the weights for plotting
    [Wa, Wi, Wo] = EcogONMF.get_clustering_weights(U,Ulb);
    if flag_renormalize
        Wa(Wa>=0) = Wa(Wa>=0) ./ max(Wa(Wa>=0),[],'all');
        Wa(Wa<=0) = Wa(Wa<=0) ./ max(abs(Wa(Wa<=0)),[],'all');
        Wi = Wi./ max(Wi, [], 'all');
        Wo = Wo./ max(Wo, [], 'all');
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

    % get areas in each RoI
    area_inds = EcogPreProcessing.get_elecs_in_area(coords.areas);

    % make a new figure;
    SCR_SZ = get(0, 'Screensize');
    SCR_SZ(end) = SCR_SZ(end)/3;
    fig = figure('Position',SCR_SZ);
    nrows = length(clusters_to_show);
    ncols = 3;
    tiledlayout(nrows,ncols, 'Padding', 'none', 'TileSpacing', 'compact');

    % set ylims for component plots
    V_min = min(V(clusters_to_show,:),[],'all');
    V_max = max(V(clusters_to_show,:),[],'all');

    % loop over cluster comps
    for i = 1:length(clusters_to_show)
        % cluster number
        k = clusters_to_show(i);
        % plot the temporal comp
        nexttile(3*i-2);
        plot(time_start_window,...
            V(k,:),...
            'LineWidth', 2);
        ylim([V_min, V_max]);
        if i==1;
            xlabel('Time (sec)');
            ylabel('Normalized Directed Connectivity');
        end
        % plot the weights
        nexttile(3*i-1);
        VT.PlotWeightedBrain(coords.(coords_to_use),...
                             Wa(:,k),...
                             'ElecRegions', coords.areas,...
                             'Sigma', 5,...
                             'Mode', 'SimpleProject',...
                             'cmap', cmap, 'clim', clim);
        % plot SI vs Inflow
        nexttile(3*i);
        Inflow_ = Wi(area_inds.(Inflow_RoI), k);
        SI_ = SI(area_inds.(Inflow_RoI));
        [r, p] = corrcoef(SI_, Inflow_);
        r = r(1,2);
        p = p(1,2);
        plot(SI_, Inflow_, 'k*', 'LineWidth', 2, 'MarkerSize',8);
        hold on;
        rline = polyfit(SI_, Inflow_,1);
        x = linspace(min(SI_)-0.01,max(SI_)+0.01,100);
        y = polyval([rline(1),rline(2)], x);
        plot(x, y, 'b--', 'LineWidth', 2);
        xlim([min(SI_)-0.05,max(SI_)+0.05])
        ylim([min(Inflow_)-0.05,max(Inflow_)+0.05])
        xlabel('suppressin index');
        ylabel('Inflow weight');
        title(sprintf('SI v Inflow, r=%1.2f, p=%1.3f', r, p));
    end
end
