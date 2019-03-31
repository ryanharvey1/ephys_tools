function str = Events2CSV(ev_fname,csv_fname)
%
% str =  Events2CSV(ev_fname,csv_fname)
%
% convert a Cheetah event file ev_fname to CSV format (for excel spreadsheet) 
% with 3 ascii columns separated by commas (csv_fname):
%
% timestamps (0.1 msec units) ,  EventCode (TTL input) , Event String (console input and TTL input)
%
% Returns a string cell array with one csv line per record.
%
% PL 2001

[evtsd, ts, param1, param2, param3, evstring] = LoadEventFlags(ev_fname);

% open output csv file
[fid, msg] = fopen(csv_fname, 'w');
if(fid == -1)
    error(msg);
end

ttl = Data(evtsd);

for ii=1:length(ts)
    str{ii} = sprintf('%12d, %d, "%s"',floor(ts(ii)), ttl(ii), evstring{ii});
    disp(str{ii});
    fprintf(fid,'%s\n',str{ii});
end

fclose(fid);



