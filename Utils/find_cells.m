function [cell_index]=find_cells(data,tetrode,cell)
% find_cells: locates cell index from ephys_tools standard
%
%   Input:
%           data: ephys_tools data stucture
%           tetrode: tetrode number (single number or vector)
%           cell: cell number (single number or vector)
%
%   Output:
%           cell_index: row number(s) location of cells
%
% Note: If tetrode is the only input, find_cells will return all cells on
% the supplied tetrode number(s)
%
% Ryan Harvey 2019

if exist('cell','var')
    for i=1:length(cell)
        cells_to_find{i} = ['TT',num2str(tetrode(i)),'.mat',num2str(cell(i))];
    end
    
    cell_list=strcat(data.spikesID.TetrodeNum,num2str(data.spikesID.CellNum));
    
    cell_index=find(ismember(strrep(cell_list,' ',''),cells_to_find))';
else
    for i=1:length(tetrode)
        cells_to_find{i} = ['TT',num2str(tetrode(i)),'.mat'];
    end
    cell_list=strcat(data.spikesID.TetrodeNum);
    cell_index=find(ismember(strrep(cell_list,' ',''),cells_to_find))';
end
end