function conditions=get_conditions
% get_conditions: list available standardized conditions
%
% you can add more conditions by editing conditions.csv in the
% preprocessing folder
%
%Ryan Harvey 2019
opts = delimitedTextImportOptions("NumVariables", 1);
opts.DataLines = [1, Inf];
opts.Delimiter = "";
opts.VariableNames = "rotatec";
opts.VariableTypes = "string";
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
conditions = readtable("conditions.csv", opts);
conditions = table2cell(conditions);
numIdx = cellfun(@(x) ~isnan(str2double(x)), conditions);
conditions(numIdx) = cellfun(@(x) {str2double(x)}, conditions(numIdx));
end

