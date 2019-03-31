function projresult=ParseAllenMouseBrainConnectivity(filename)
%% Import data from text file.

%% Initialize variables.
delimiter = ',';
startRow = 2;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,3,4,8,12,13,17]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,3,4,8,12,13,17]);
rawStringColumns = string(raw(:, [2,5,6,7,9,10,11,14,15,16]));


%% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
for catIdx = [1,2,3,6,7,10]
    idx = (rawStringColumns(:, catIdx) == "<undefined>");
    rawStringColumns(idx, catIdx) = "";
end

%% Create output variable
projresult = table;
projresult.id = cell2mat(rawNumericColumns(:, 1));
projresult.transgenicline = categorical(rawStringColumns(:, 1));
projresult.productid = cell2mat(rawNumericColumns(:, 2));
projresult.structureid = cell2mat(rawNumericColumns(:, 3));
projresult.structureabbrev = categorical(rawStringColumns(:, 2));
projresult.structurename = categorical(rawStringColumns(:, 3));
projresult.name = rawStringColumns(:, 4);
projresult.injectionvolume = cell2mat(rawNumericColumns(:, 4));
projresult.injectionstructures = rawStringColumns(:, 5);
projresult.gender = categorical(rawStringColumns(:, 6));
projresult.strain = categorical(rawStringColumns(:, 7));
projresult.sum = cell2mat(rawNumericColumns(:, 5));
projresult.structurecolor = cell2mat(rawNumericColumns(:, 6));
projresult.numvoxels = rawStringColumns(:, 8);
projresult.injectioncoordinates = rawStringColumns(:, 9);
projresult.selected = categorical(rawStringColumns(:, 10));
projresult.experiment_page_url = cell2mat(rawNumericColumns(:, 7));

%% create simple to read stats

regions=cellstr(unique(projresult.structureabbrev));
site=[];
for i=1:length(regions)
    idx=find(contains(cellstr(projresult.structureabbrev),regions{i}));
    
    cell2mat(strfind(cellstr(projresult.structureabbrev),regions{1}))
    
    index = cellfun(@(x) x==1, strfind(cellstr(projresult.structureabbrev),regions{1}), 'UniformOutput', 0)


    for ii=idx'
        site=[site;extractBetween(projresult.injectionstructures(ii),'"abbreviation"=>"','",')];
    end
    
    a=unique(site,'stable')
b=cellfun(@(x) sum(ismember(site,x)),a,'un',0)
    
end








end
