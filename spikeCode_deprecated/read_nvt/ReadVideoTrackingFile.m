% Matlab Program to read a Neurolynx Video tracking Data File
% Written by Susan Schwarz, Academic Computing dartmouth Collelge 9/2008

HEADER_SIZE=16384;
RECORD_SIZE=1828;
NLX_VTREC_NUM_POINTS=400;
NLX_VTREC_NUM_TARGETS=50;
[name,filepath]=uigetfile('.nvt','Select a Neuralynx video tracking file');
filename=fullfile(filepath,name);
finfo=dir(filename);
filesize=finfo.bytes;
nrecs= floor((filesize-HEADER_SIZE)/RECORD_SIZE);
fid =fopen(filename,'r','a');
header=fread(fid,16384,'uchar=>uchar');
swstx=zeros(1,nrecs,'uint16');
swid=zeros(1,nrecs,'uint16');
swdata_size=zeros(1,nrecs,'uint16');
qwTimeStamp=zeros(1,nrecs,'uint64');
%dwPoints=zeros(400,nrecs,'uint32');
dwPoints=zeros(400,'uint32');
sncrc=zeros(1,nrecs,'int16');
dnextracted_x=zeros(1,nrecs,'int32');
dnextracted_y=zeros(1,nrecs,'int32');
dnextracted_angle=zeros(1,nrecs,'int32');
dntargets=zeros(50,nrecs,'int32');
targets=zeros(1,50);

% read in all of the data
for i=1:nrecs
   if (mod(i,5000) ==0)
      step=sprintf('reading record # %d of %d \n',i,nrecs);
      disp(step);
      datestr(now);
   end
   swstx(1,i)=fread(fid,1,'*uint16');
   swid(1,i)=fread(fid,1,'*uint16');
   swdata_size(1,i)=fread(fid,1,'*uint16');
   qwTimeStamp(1,i)=fread(fid,1,'*uint64');
   %dwPoints(i,1:NLX_VTREC_NUM_POINTS)=fread(fid,NLX_VTREC_NUM_POINTS,'*uint32');
   dwPoints(1:NLX_VTREC_NUM_POINTS)=fread(fid,NLX_VTREC_NUM_POINTS,'*uint32');
   sncrc(1,i)=fread(fid,1,'*int16');
   dnextracted_x(1,i)=fread(fid,1,'*int32');
   dnextracted_y(1,i)=fread(fid,1,'*int32');
   dnextracted_angle(1,i)=fread(fid,1,'*int32');
   targets= fread(fid,NLX_VTREC_NUM_TARGETS,'*int32');
   dntargets(:,i)=targets;
end
fclose(fid);
disp('finished reading video tracking file\n');
% Extract the x,y positions for red and green from the target data
% open the file for the x,y target data
loc=strfind(filename,'.nvt');
target_data_filename=strcat(filename(1:loc-1),'.txt');
fid=fopen(target_data_filename,'w');
disp('starting to process video tracking data file.');
for i =1:nrecs
   if (mod(i,5000) ==0)
      step=sprintf('processing record # %d of %d \n',i,nrecs);
      disp(step);
      datestr(now);
   end
   red_x=0;
   red_y=0;
   green_x=0;
   green_y=0;
   for j=1:4  % only need to use the first 4 values
      [red, green, blue, intensity, x, y] = vt_target_decode(dntargets(j,i)) ;
      %output=sprintf('red=%d,green=%d,blue=%d,intensity=%d,x=%d,y=%d',red, green, blue, intensity, x, y);
      %disp(output);
      if (red ==1) 
	 red_x=x;
	 red_y=y;
      elseif (green==1)
	  green_x=x;
	  green_y=y;
      end
   end
   fprintf(fid,'%u\t%d\t%d\t%d\t%d\n',qwTimeStamp(1,i),red_x,red_y,green_x,green_y);
end
fclose(fid);
str=sprintf('Finished writing file: %s \n', target_data_filename);
disp(str);




