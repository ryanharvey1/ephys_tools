% Process DLC for postprocess

function [ts, x, y, angles] = process_DLC_for_ephys(varargin)
% calculates behavior correlates for data tracked using DLC. Accomodates
% two behavioral pipelines from BClark lab (custom Open Ephys and Neuralynx). 
% 
% Assumes folder/file structure for ephys_tools. Requires animal metadata
% folder with led positions (either front,back,left,right relative to animal head).
% DLC tracking points should match up with the strings indicated in
% metadata (i.e. AnimalMetadata.ExtracellEphys.led_position.front = 'red'
% and in 'red*' in DLC header). 
% 
% Laura Berkowitz 2020
% updated with Neuralynx tracking LB March 2021

% Inputs 
p = inputParser;
p.addParameter('basepath',pwd);
p.addParameter('metadata','AnimalMetaData');
p.parse(varargin{:});

basepath = p.Results.basepath;
metadata = p.Results.metadata;

% Load animal metadata
pathparts = strsplit(basepath,filesep);
metadata_path = [pathparts{1},filesep,pathparts{2},filesep,metadata];
load([metadata_path,filesep,pathparts{4},'_metadata.mat']);

% grab led position
front = AnimalMetadata.ExtracellEphys.led_position.front;
back = AnimalMetadata.ExtracellEphys.led_position.back;
left = AnimalMetadata.ExtracellEphys.led_position.left;
right = AnimalMetadata.ExtracellEphys.led_position.right;

% Assign led color to location
if sum(ismember([front{1}],'nan')) == 3
    % tracking points are right/left
    bodypart1 = left{1};   
    bodypart2 = right{1};
else 
    % tracking points are back/front
    bodypart1 = front{1};   
    bodypart2 = back{1};
end


% Find DLC excel
file = table2cell(struct2table(dir([basepath,filesep,'*.csv'])));
dlc_idx = contains(file(:,1),'DLC') & ~contains(file(:,1),'above95th');

if sum(dlc_idx) == 0
    disp('No DLC tracking found')
    return
end

dlc_file = file(dlc_idx,:);

% Pull XY for red/yellow
[header,tsxy]= open_dlc(dlc_file);

% get x/y coordinates from bodyparts
bodypart1 = tsxy(:,contains(header(1,:),bodypart1) & ~contains(header(2,:),'likelihood'));
bodypart2 = tsxy(:,contains(header(1,:),bodypart2) & ~contains(header(2,:),'likelihood'));

% median between x and y is center of the head, so lets make those our new
% x/y values. 
x = nanmedian([bodypart1(:,1) bodypart2(:,1)],2)';
y = nanmedian([bodypart1(:,2) bodypart2(:,2)],2)';

 % Compute head angle
angles = rad2deg(XYangleLED(bodypart1(:,1),bodypart1(:,2),bodypart2(:,1),bodypart2(:,2)))';
if sum(ismember([front{1}],'nan')) == 3
    angles = wrapTo360(angles-90); % account for markers on side of head instead of front/back
end

% Sync start of spikes with video 
if ~isempty(dir([basepath,filesep,'*record_ts.csv']))
    % get file directory info
    ts_file = dir([basepath,filesep,'*record_ts.csv']);
    vid_file = dir([basepath,filesep,'*video_ts.csv']);
    % load record/start file
    rec_ts = readtable(ts_file.name);
    % load video timestamps and remove the first row (it's initialized to
    % zero before recording starts)
    vid_ts = readtable(vid_file.name); vid_ts(1,:) = []; 
    % find difference between recording start and video timestamp start
    vid_offset = vid_ts.video_ts(1,1) - rec_ts.start_record;
    % add video offset to ts. 
    ts = (vid_ts.video_ts - vid_ts.video_ts(1,1))+ vid_offset; 
    % convert to microseconds 
    ts = [ts*10^6]';
    
% for neuralynx bitfields     
elseif ~isempty(dir([basepath,filesep,'*.ntt']))
    ts = Nlx2MatVT([basepath,filesep,'VT1.nvt'],[1,0,0,0,0,0],0,1);
    % remove trailing frames (Neuralynx capture has delay in closing video
    % capture, so frames may bigger). 
    ts = ts(1:length(angles));
else 
    % For now, just base timestamps off of DLC frames. These sessions are noted and will likely not be used given the ambiguous offset.  
    ts = [0:1/30:tsxy(end,1)/30]*10^6;
end

end

function [header,tsxy]= open_dlc(dlc_file)
% dlc_file is cell array of 'dir' output.
% for example,
% file = table2cell(struct2table(dir([path,filesep,'*.csv'])));
% where path is the path to the dlc output

% load header
fileID = fopen(fullfile(dlc_file{2},dlc_file{1}),'r');
dataArray = textscan(fileID, '%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]',...
    3-2+1, 'Delimiter',',', 'TextType', 'string', 'HeaderLines',...
    2-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
header = [dataArray{1:end-1}];
clearvars fileID dataArray ans;

% load data
fileID = fopen(fullfile(dlc_file{2},dlc_file{1}),'r');
dataArray = textscan(fileID, '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]',...
    'Delimiter', ',','TextType', 'string', 'EmptyValue', NaN,...
    'HeaderLines' ,4-1,'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
tsxy = [dataArray{1:end-1}];

% Turn low liklihood point into NAN and smooth remaining coords
idx = tsxy(:,contains(header(2,:),'likelihood'))<.95;
xloc = find(contains(header(2,:),'x'));
yloc = find(contains(header(2,:),'y'));
for l = 1:size(idx,2)
    tsxy(idx(:,l),xloc(l):yloc(l))=NaN;

end

end