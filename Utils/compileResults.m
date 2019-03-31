function data=compileResults(ProcessedData) 
% compileResults: compiles measures and cell ids into a stucture. This function also
% displays helpful stats: n animals, animal ids, n recording sessions,
% n clusters per animal.
% 
%   Input: ProcessedData: path to processed data .mat files. These .mat
%   files should store the standard bclarklab data structure
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

%% EXTRACT MEASURES AND IDS
data=[];
for i=1:length(sessions)
   load(sessions{i},'measures','rat','spikesID') 
   [r,c,d]=size(measures);
   tempfill=nan(r,c,nsessions);
   tempfill(:,:,1:d)=measures;
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

load(sessions{1},'varnames')
data.varnames=varnames;

end