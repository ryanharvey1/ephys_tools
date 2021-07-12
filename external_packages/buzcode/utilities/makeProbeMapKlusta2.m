function [] = makeProbeMapKlusta2(folder, xmlfile)
%% this function takes an xml file with spike groups already assigned and 
%% generates *.prb files for each shank for spike detection/masking/clustering

%% IMPORTANT: this function ASSUMES that the order of the channels for each shank
%% in the spike groups map onto a specific geometrical arrangement of channels on the shank
%% starting from the top left recording site, and moving right then downward to make the 
%% nearest neighbor graph
%make folders
%make proper prb in each folder
%put prm in each folder

% parameters = LoadParameters([folder '/' xmlfile]);
% warning off
% 
% for shank = 1:length(parameters.spikeGroups.groups)

if ~exist('folder','var')
    folder = cd;
end
if ~exist('xmlfile','var')
%     d = dir([folder,'/*.xml_forklusta']);%this doesn't read in right... 
%     if isempty(d)
        d = dir([folder,'/*.xml']);
%     end
    xmlfile = d(1).name;
end

parameters = LoadPar([folder '/' xmlfile]);
warning off
if exist(fullfile(folder,'bad_channels.txt'),'file')
    badchannels = ReadBadChannels(folder);
else
    badchannels = [];
end
spkgroupnums = matchSpkGroupsToAnatGroups(parameters);


for shi = 1:length(parameters.SpkGrps)
    
    shank = spkgroupnums(shi);%the orignal name of the shank, not just sequential number
    % make a folder for each directory
    if 1
%     if exist([folder '/' num2str(shank)],'dir')
%        disp([folder '/' num2str(shank) ' already exists, not over-writing'])
%     else


    %     channels = parameters.spikeGroups.groups{shank};
        channels = parameters.SpkGrps(shi).Channels;
        
        if ~isempty(setdiff(channels,badchannels))
            mkdir([folder '/' num2str(shank)]);

            c=1;
            for i=1:length(channels)-1
                for j=1
                l(c,:) = [channels(i),channels(i+j)];
                c=1+c;
                end
            end    
            for i=1:length(channels)-2
                for j=2 
                l(c,:) = [channels(i),channels(i+j)];
                c=1+c;
                end
            end
            for i=1:length(channels)-3
                for j=3
                l(c,:) = [channels(i),channels(i+j)];
                c=1+c;
                end
            end
            list = l;    

            badrows = [];
            for a = 1:length(badchannels)
                [t,~] = find(list==badchannels(a));
                badrows = cat(1,badrows,t);
            end
            list(badrows,:) = [];
            if ~isempty(badrows)
                1;
            end

            s=['channel_groups = {\n' num2str(shank) ': {\n'];

            s=[s, '\t"channels": [\n' ];
            for i =1:length(channels)-1
                if ~ismember(channels(i),badchannels)%don't record bad channels
                    s=[s, '' num2str(channels(i)) ', ' ];
                end
            end
            s=[s, '' num2str(channels(i+1))];
            s=[s, '],\n' ];

        %     !handle badchannels here!... 
            s=[s, '\t"graph": [\n' ];
            for i =1:length(list)-1
                s=[s, '\t(' num2str(list(i,1)) ', ' num2str(list(i,2)) '),\n'];
            end
            s=[s, '\t(' num2str(list(i+1,1)) ', ' num2str(list(i+1,2)) ')\n'];
            s=[s, '\t],\n' ];

            s=[s, '\t"geometry": {\n' ];
            for i =1:length(channels)
                if ~ismember(channels(i),badchannels)
                    x = length(channels)-i;
                    y = -i*1;
                    if mod(i,2)
                        x = -x;
                    end
                    s=[s, '\t' num2str(channels(i)) ': [' num2str(x) ', ' num2str(y) '], \n'];
                end
            end

    %         for i =1:length(channels)
    %             if ~ismember(channels(i),badchannels)
    %                 x = pn*.5*length(channels)-i;
    %                 y = -i*1;
    %                 if mod(i,2)
    % %                    pn = -1;
    %                    s=[s, '\t' num2str(channels(i)) ': [' num2str(-x) ', ' num2str(y) '], \n'];
    %                 else
    % %                    pn = 1;
    %                    s=[s, '\t' num2str(channels(i)) ': [' num2str(x) ', ' num2str(y) '], \n'];
    %              end
    % %                 s=[s, '\t' num2str(channels(i)) ': [' num2str(pn*(-(2+1-i)*2)) ', ' num2str(-i*1) '], \n'];
    %             end
    %         end
            s=[s, '\t},\n},\n}\n' ];

            % write .prb files to each shank folder
            fid = fopen([folder '/' num2str(shank) '/' num2str(shank) '.prb'],'wt');
            fprintf(fid,s);
            fclose(fid); 

            % write .prm file as well
            fid = fopen('D:/Users/BClarkLab/ephys_tools/preprocessing/template.prm');
    %         fid = fopen('/mnt/packrat/userdirs/brendon/ProcessingFiles/TemplateOrig.prm');
            i = 1;
            tline = fgetl(fid);
            generic{i} = tline;
            while ischar(tline)
                i = i+1;
                tline = fgetl(fid);
                generic{i} = tline;
            end
            fclose(fid);
        %         'filebase = ''' folder '/' num2str(shank) '/' xmlfile(1:end-4) '_sh' num2str(shank) '''\n',...
            ss = ['experiment_name = ''' xmlfile(1:end-4) '_sh' num2str(shank) '''\n',...
                'filebase = experiment_name\n',...
                'prb_file = ''' num2str(shank) '.prb''\n',...
                '\n',...
                'traces = dict(\n',...
                '    raw_data_files=[''' folder  '/' xmlfile(1:end-4) '.dat'  '''],\n',...
                '    n_channels = ' num2str(parameters.nChannels) ',\n',...            
                ];

            fid = fopen([folder '/' num2str(shank) '/' xmlfile(1:end-4) '_sh' num2str(shank) '.prm'],'wt');
            fprintf(fid,ss);
            for i = 1:numel(generic)
                if generic{i+1} == -1
                    fprintf(fid,'%s', generic{i});
                    break
                else
                    fprintf(fid,'%s\n', generic{i});
                end
            end
            clear s l list c ss
        end
    end
end
disp('Shank directories and probe maps made')



function spkgroupnums = matchSpkGroupsToAnatGroups(par)
% find the original number of the shank, even if a shank was removed -
% using matches to anat group numbers.  Without this elec groups are just
% named in sequence regardless of skips of actual shanks.  Assumes 1:1
% matching of channels to groups

nanat = length(par.AnatGrps);
nspk = par.nElecGps;
spkgroupnums = zeros(nspk,1);
anatfirsts = zeros(nanat,1); %make a list of the first chan of each anat group

for a = 1:nspk
   spkchans = par.SpkGrps(a).Channels;
   for b = 1:nanat
       lia = ismember(par.AnatGrps(b).Channels,par.SpkGrps(a).Channels);
       if any(lia)
%            lia;
           spkgroupnums(a) = b;
           continue
       end
   end
end
1;
% 
% for a = 1:nanat;
%    anatfirsts(a) = par.AnatGrps(a).Channels(1);
% end
% 
% for a = 1:nspk
%    spkchans = par.SpkGrps(a).Channels;
%    spkgroupnums(a) = find(ismember(anatfirsts,spkchans));
% end
% 