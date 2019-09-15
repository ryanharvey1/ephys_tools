function data=compileResults(ProcessedData) 
% compileResults: compiles measures and cell ids into a stucture. This function also
% displays helpful stats: n animals, animal ids, n recording sessions,
% n clusters per animal.
% 
%   Input: ProcessedData: path to processed data .mat files. These .mat
%   files should store the standard ephys_tools data structure
%
%   Output: data: structure of measures organized by rat id
%
% Ryan Harvey
%
cd(ProcessedData)
sessions=dir('*.mat');
sessions={sessions.name}';

%% READ OUT BASIC INFO ABOUT DATA SET
tempsess=split(sessions,'_');
rats=unique(tempsess(:,1));
% n animals
disp([num2str(length(rats)),' animals'])
% id
disp(rats')

%% GET NUMBER OF RECORDING SESSIONS, MAX NUMBER OF SESSIONS PER RECORDING SESSION,
% AND DELETE THE FILE NAMES OF SESSIONS THAT WERE TOO SHORT. 
warning off
sessiontodelete=[];
for i=1:length(sessions)
   load(sessions{i},'measures')
   if exist('measures','var')
        nsessions(i,1)=size(measures,3);
   else
       sessiontodelete=[sessiontodelete;i];
   end
   clear measures
end
warning on
sessions(sessiontodelete)=[];
disp([num2str(length(nsessions)),' Recording sessions'])
nsessions=max(nsessions);

%% get varnames (for cases where variables were recently added)
for i=1:length(sessions)
    load(sessions{i},'varnames')
    all_vars{i}=varnames;
end
all_varnames=unique([all_vars{:}],'stable');
c=length(all_varnames);

%% EXTRACT MEASURES AND IDS
data=[];
for i=1:length(sessions)
   load(sessions{i},'measures','rat','spikesID','varnames') 
   [r,~,d]=size(measures);
   tempfill=nan(r,c,nsessions);
   % check and pull out each variable 
   for v=1:c
       if contains(all_varnames{v},varnames)
           tempfill(:,v,1:d)=measures(:,ismember(varnames,all_varnames{v}),1:d);
       else
           continue
       end
   end
   if ~isfield(data,rat)
       data.(rat).measures=tempfill;
       data.(rat).id=[cellstr(repmat(sessions{i},r,1)),spikesID.TetrodeNum,cellstr(num2str(spikesID.CellNum))];
   else
       data.(rat).id=[data.(rat).id;[cellstr(repmat(sessions{i},r,1)),spikesID.TetrodeNum,cellstr(num2str(spikesID.CellNum))]];
       data.(rat).measures=cat(1,data.(rat).measures,tempfill);
   end
   clear measures rat tempid tempfill
end

%% DISP CLUSTER NUMBER PER RAT AND OVERALL 
for i=1:length(rats)
    nclust(i,1)=size(data.(rats{i}).measures,1);
    disp([rats{i},': ',num2str(nclust(i,1)),' Clusters'])
end
disp([num2str(sum(nclust)),' Clusters Overall'])

data.varnames=all_varnames;

end