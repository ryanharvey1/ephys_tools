function Sep_ntt_by_event(raw_dir)


d = dir([raw_dir '\*.ntt']);
for i = 1:length(d)
    disp(['Filtering: ' raw_dir '\' d(i).name]);
    [path, name, ext] = fileparts(d(i).name);
    temp = sscanf(name, '%c%c%f');
    TTnum = temp(3);
    outfile = [raw_dir filesep 'TT' num2str(TTnum) '.ntt'];
    disp('FILTERING BY EVENT');
    event_filter(outfile, outfile);
    
end
end
function event_filter(tetrode_file,outfile)

[path, ~, ~] = fileparts(tetrode_file);
[ timestamps,StartofRec,EndofRec ] = EventSplit( path );


[num_recs_all, ~, ~] = GetNlxFileProperties_v4(tetrode_file);
% disp(['File contains: ' num2str(num_recs_all) ' spikes']);


if num_recs_all<2000000
    % KEVIN'S FILE LOADER FOR SMALLER FILES
    % set up the fields used by Nlx2MatSpike
    FieldSelection(1) = 1;
    FieldSelection(2) = 1;
    FieldSelection(3) = 1;
    FieldSelection(4) = 1;
    FieldSelection(5) = 1;
    ExtractHeader = 1; % yes, do
    ExtractionMode = 1;  % Extract Record Index Range
    
    % note we subtract one from the look-up indices because Nlx2MatSpike
    % uses C-style indices (starting from 0)
    [ts, ScNumbers, CellNumbers, Params, wv_data, NLX_header] = ...
        Nlx2MatSpike_v4(tetrode_file, FieldSelection,ExtractHeader, ExtractionMode);
else
    disp('Using LoadTT_NeuralynxNT to load files because file is very large');
    [ts, wv_data] = LoadTT_NeuralynxNT(tetrode_file,[0 num_recs_all], 4);
    ts = ts'*100;  % reformat and convert to microseconds
    wv_data = permute(wv_data, [3 2 1]);  % reformat
    
    % now get header
    FieldSelection(1) = 0;
    FieldSelection(2) = 0;
    FieldSelection(3) = 0;
    FieldSelection(4) = 0;
    FieldSelection(5) = 0;
    ExtractHeader = 1; % yes, do
    ExtractionMode = 2;  % Extract Record Index Range
    

    NLX_header = ...
        Nlx2MatSpike_v4(tetrode_file, FieldSelection,ExtractHeader, ExtractionMode,[0 1]);
    %         end
end

tstemp=[];
% ScNumberstemp=[];
% CellNumberstemp=[];
% Paramstemp=[];
wv_datatemp=[];
for event=1:length(StartofRec)
    tstemp=[tstemp,ts(:,ts>StartofRec(event) & ts<EndofRec(event))];
%     ScNumberstemp=[ScNumberstemp,ScNumbers(:,ts>StartofRec(event) & ts<EndofRec(event))];
%     CellNumberstemp=[CellNumberstemp,CellNumbers(:,ts>StartofRec(event) & ts<EndofRec(event))];
%     Paramstemp=[Paramstemp,Params(:,ts>StartofRec(event) & ts<EndofRec(event))];
    wv_datatemp=cat(3,wv_datatemp,wv_data(:,:,ts>StartofRec(event) & ts<EndofRec(event)));

end

% 	if write_file_flag
num_recs = length(tstemp);
% clear wv_data;  % temporary fix for memory management issue with Nlx2MatSpike
scn = zeros(1,num_recs);
pars = zeros(8,num_recs);
cn = zeros(1,num_recs);

% WRITE RECORDS
% set up the fields used by Mat2NlxTT_v4
wFieldSelection(1) = 1;  % timestamps
wFieldSelection(2) = 1;  % Sc Numbers
wFieldSelection(3) = 1;  % Cell Numbers
wFieldSelection(4) = 1;  % Params
wFieldSelection(5) = 1;  % Data Points
wFieldSelection(6) = 1;  % Header
AppendFile = 0; % no, don't append
ExtractMode = 1;  % Write all records

Mat2NlxTT_v4(outfile, AppendFile, ExtractMode, 1, num_recs, wFieldSelection, tstemp, scn, cn, pars, wv_datatemp, NLX_header);
end




