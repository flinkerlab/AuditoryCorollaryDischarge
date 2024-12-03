function [colors] = get_per_task_colors()
    % returns a struct that is a lookup for colors for each task
    colors = struct();
    colors.PicN    = [0.2627    0.3608    0.9804];
    colors.VisRead = [0.3373    0.7020    0.7294];
    colors.AudRep  = [0.3922    0.8314    0.0745];
    colors.AudN    = [0.8706    0.1373    0.5882];
    colors.SenComp = [0.6350    0.0780    0.1840];
    colors.AudRep_passive = [0.9290 0.6940 0.1250];
end
