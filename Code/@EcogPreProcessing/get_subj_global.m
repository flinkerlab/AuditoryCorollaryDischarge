function [glob] = get_subj_global(obj)
	glob.subj = obj.subj;
	glob.rootdir = obj.root_path;
    if strcmp(obj.task, 'AudRep_passive')
        glob.task = 'AudRep_passive';
        % glob.task = 'AudRepPassive';
    else
        glob.task = obj.task;
    end
    glob.SJdir = [obj.root_path, filesep, 'Data', filesep, obj.subj];
	glob.ANdir = [glob.SJdir, filesep, 'analysis', filesep, obj.task];
	glob.DTdir = [glob.SJdir, filesep, 'data', filesep, obj.task];
	glob.Eventsdir = [glob.ANdir, filesep, 'Events.mat'];
	glob.CSVdir = [glob.SJdir, filesep, 'analysis', filesep,...
				   'coordinates.csv'];
	subj_glob_file = [glob.ANdir, filesep, 'subj_globals.mat'];
	if exist(subj_glob_file, 'file')
		glob.bad_elecs = load(subj_glob_file, 'bad_elecs').bad_elecs;
		glob.srate = load(subj_glob_file, 'srate').srate;
		glob.ANsrate = load(subj_glob_file, 'ANsrate').ANsrate;
	end
    active_elecs_file = [glob.SJdir, filesep, 'analysis',...
                         filesep, 'Active_elecs.mat'];
    glob.selected_elecs = obj.selected_elecs;
    if obj.flag_load_only_active_elecs
        if exist(active_elecs_file, 'file')
            SE = load(active_elecs_file);
            glob.selected_elecs.active_elecs = SE.active_elecs;
            glob.selected_elecs.active_elecs_labels = SE.active_elecs_labels;
        else
            glob.selected_elecs = [];
        end
    end
end

