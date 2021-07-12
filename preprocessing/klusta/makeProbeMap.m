function [] = makeProbeMap(folder,varargin)
%% this function takes an xml file with spike groups already assigned and
%% generates *.prb files for each shank for spike detection/masking/clustering

%% IMPORTANT: this function ASSUMES that the order of the channels for each shank
%% in the spike groups map onto a specific geometrical arrangement of channels on the shank
%% starting from the top left recording site, and moving right then downward to make the
%% nearest neighbor graph

% Pulls parameters from the xml. Can add xcoords,ycoords for spatial
% position of probe. Default coords are compatible with cambridge neurotech assy156.

% Adapted for ephys_tools from David Tingley's open_ephys_pipeline github repo
p = inputParser;
p.addParameter('xcoords',[0,70,5,65,...% shank A
    10,60,15,55,...
    20,50,25,45,...
    30,40,35,35,...
    0+250,70+250,5+250,65+250,... % shank B
    10+250,60+250,15+250,55+250,...
    20+250,50+250,25+250,45+250,...
    30+250,40+250,35+250,35+250,...
    0+500,70+500,5+500,65+500,... % shank C
    10+500,60+500,15+500,55+500,...
    20+500,50+500,25+500,45+500,...
    30+500,40+500,35+500,35+500,...
    0+750,70+750,5+750,65+750,... % shank D
    10+750,60+750,15+750,55+750,...
    20+750,50+750,25+750,45+750,...
    30+750,40+750,35+750,35+750]');
p.addParameter('ycoords', [0,-25,-40,-65,...% shank A
    -80, -105, -120, -145,...
    -160,-185,-200,-225,...
    -240,-265,-280,-305,...
    0,-25,-40,-65,... % shank B
    -80, -105, -120, -145,...
    -160,-185,-200,-225,...
    -240,-265,-280,-305,...
    0,-25,-40,-65,... % shank C
    -80, -105, -120, -145,...
    -160,-185,-200,-225,...
    -240,-265,-280,-305,...
    0,-25,-40,-65,... % shank D
    -80, -105, -120, -145,...
    -160,-185,-200,-225,...
    -240,-265,-280,-305]');
p.addParameter('dat_file','filtered_CAR.dat')
p.parse(varargin{:});

% Load xml from file.
d   = dir([folder,filesep,'*.xml']);
parameters = LoadXml(fullfile(folder,d(1).name));
[~,basename] = fileparts(folder);

% Unpack xml and probe geometry
xcoords = p.Results.xcoords;
ycoords = p.Results.ycoords;
Nshanks = parameters.nElecGps;
fs = parameters.SampleRate;
Nchannels = length(xcoords);
kcoords = (reshape(repmat(1:Nshanks, Nchannels/Nshanks, 1), Nchannels, 1));
DTYPE = 'int16';
dat_file = p.Results.dat_file;

% Loop through shanks
for shank = 1:parameters.nElecGps
    bad_channels = parameters.AnatGrps(shank).Skip';
    channels = parameters.AnatGrps(shank).Channels;
    positions = [xcoords(kcoords == shank) ycoords(kcoords == shank)];
    
    try
        % make a folder for each directory
        if ~exist([folder,filesep,'klusta_' ,num2str(shank)])
            disp(['working on spike group #' num2str(shank)])
            mkdir([folder,filesep,'klusta_' , num2str(shank)]);
            
            channels = parameters.AnatGrps(shank).Channels;
            
            % make adjacency graph.
            c=1;
            for i=1:length(channels)-1
                for j=1
                    l(c,:) = [channels(i),channels(i+j)];
                    c=1+c;
                end
            end
            
            if length(channels) < 32
                
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
            end
            list = l;
            
            % write first line of text file
            s=['channel_groups = {\n' num2str(shank) ': {\n'];
            
            % second line
            s=[s, '\t"channels": [\n' ];
            
            % loop through channels to get channel list
            for i =1:length(channels)-1
                s=[s, '' num2str(channels(i)) ', ' ];
            end
            
            % close list
            s=[s, '' num2str(channels(i+1))];
            s=[s, '],\n' ];
            
            % start writing graph section
            s=[s, '\t"graph": [\n' ];
            for i =1:length(list)-1
                s=[s, '\t(' num2str(list(i,1)) ', ' num2str(list(i,2)) '),\n'];
            end
            %     s=[s, '\t(' num2str(list(i+1,1)) ', ' num2str(list(i+1,2)) ')\n'];
            s=[s, '\t],\n' ];
            
            
            % write geometry section
            s=[s, '\t"geometry": {\n' ];
            for i =1:length(channels)
                s=[s, '\t' num2str(channels(i)) ': [' num2str(positions(i,1)) ', ' num2str(positions(i,2)) '], \n'];
            end
            s=[s, '\t},\n},\n}' ];
            
            % write .prb files to each shank folder
            fid = fopen([folder,filesep,'klusta_' ,num2str(shank),filesep,num2str(shank), '.prb'],'wt');
            fprintf(fid,s);
            fclose(fid);

            %% GET TEMPLATE INFO
            % Get template for klusta parameters file (PRM)
            kp_template = import_klusta_pars;
            % Write PRM
            exp_folder = strrep(folder,'\','/');
            exp_folder_new = [exp_folder '/' 'klusta_' num2str(shank) '/' basename '_sh' num2str(shank)];
            pdir = fullfile([folder,filesep,'klusta_' , num2str(shank)]);
            exp = [basename '_sh' num2str(shank)];
            fname = fullfile(pdir,[exp '.prm']);
            fid = fopen(fname,'w');
            prb_name = [num2str(shank) '.prb'];
            for iLine = 1:numel(kp_template)
                switch iLine
                    case 1
                        fwrite(fid,sprintf([char(kp_template(iLine)) '\n'],exp_folder_new));
                    case 2
                        fwrite(fid,sprintf([char(kp_template(iLine)) '\n'],prb_name));
                    case 4
                        fwrite(fid,sprintf([char(kp_template(iLine)) '\n'],[exp_folder '/']));
                    case 5
                        fwrite(fid,sprintf([char(kp_template(iLine)) '\n'],10));
                    case 6
                        fwrite(fid,sprintf([char(kp_template(iLine)) '\n'],fs));
                    case 7
                        fwrite(fid,sprintf([char(kp_template(iLine)) '\n'],Nchannels));
                    case 8
                        fwrite(fid,sprintf([char(kp_template(iLine)) '\n'],DTYPE));
                    otherwise
                        fwrite(fid,sprintf([char(kp_template(iLine)) '\n']));
                end
            end
            fclose(fid);
            
            
        end
        
    catch
        warning(['did not work correctly for spike group ' shank]);
    end
end
end