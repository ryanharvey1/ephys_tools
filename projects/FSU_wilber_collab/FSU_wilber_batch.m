% FSU_wilber_batch

path='D:\Projects\FSU_wilber_collab\data';
cd(path)
sessions=dir;
sessions={sessions.name};
sessions(contains(sessions,[".",".."]))=[];

sessions=strcat(path,filesep,sessions);

for i=1:length(sessions)
   FSU_wilber_collab_v2(sessions{i}) 
end