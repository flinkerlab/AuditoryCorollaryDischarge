function out = Qevents(Database, ev, fields,str_fld,str_val)
	flds = fieldnames(Database); ev_sz = length(ev);
	if (nargin<5)
		str_val = ''; str_fld = 'event';
	end
	if (~isfield(Database,fields))
		error('fields must be members of Database !');
	end
	if (~isfield(Database,str_fld))
		error('str_fld must be members of Database !');
	end
	if (strcmp(flds(1),'event') == 0)
		error('First field is not ''event'' !');
	end
	if (~iscell(fields))
		error('fields has to be a cell array {} of the desired structure fields');
	end
	cnt = 1;
	for i=1:length(Database) % search Databse for events
		i_ev_sz = length(Database(i).event);
		% Database events starts with &ev
		if(strcmp(Database(i).event(1:min(ev_sz,i_ev_sz)), ev))
			if (~strcmp(str_fld,'event'))
				if (~strcmp(str_val,eval(['Database(i).' str_fld ';'])))
					continue;
				end
			end
			for f=1:length(fields)
				try
					eval(['out(cnt,f) = Database(i).' fields{f} ';']);
				catch ME
					out(cnt,f) = NaN;
				end
			end
			cnt = cnt+1;
		end
	end
end


