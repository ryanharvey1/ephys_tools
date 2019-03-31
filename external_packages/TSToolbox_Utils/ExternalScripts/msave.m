% msave - save an integer matrix in an ascii file
%
% msave(filename,matrix)
%

function msave(filename,matrix)

[m,n] = size(matrix);
formatstring = '%d';
for ii=2:n,
  formatstring = [formatstring,'\t%d'];
end
formatstring = [formatstring,'\n'];

outputfile = fopen(filename,'w');
fprintf(outputfile,formatstring,matrix');
fclose(outputfile);


