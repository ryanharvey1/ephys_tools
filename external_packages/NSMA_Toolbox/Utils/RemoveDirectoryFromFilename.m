function  list2 = RemoveDirectoryFromFilename(list1)
%
% Removes directories from a cell array of  filename strings
% input
%   list1 ... cell array of filenames in format e.g. D:\dir1\dir2\filename.ext
% output
%   list2 ... cell array of filenames without direcories  e.g. filename.ext
%

list2 = cell(size(list1));
for ii = 1:length(list1)
   [dir,fname,ext] = fileparts(list1{ii});
   list2{ii} = [ fname ext];
end
