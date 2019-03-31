% function aligntsp2(file, [front, back])
%
% Making whl file for two colors. Originally from 'aligntsp.m' by Tony Berenyi.
%
% tsp(t,:) = [timestamp, color1-X, color1-Y, color2-X, color2-Y, color3-X, color3-Y]
% 'front' and 'back' assigns which color (1, 2 or 3) should be used. if -1, then all rows are -1.

function aligntsp2mod(file, colIx)

if colIx(1)>0; colorf = 2*colIx(1)+[0, 1]; else colorf = [8, 8]; end
if colIx(2)>0; colorb = 2*colIx(2)+[0, 1];  else colorb = [8, 8]; end
mycolor = [colorf, colorb];

% get timestamps and positions from tsp file
tspdata=load([file '.tsp']);


%get start and end timestamp of dat file
fid=fopen([file '.meta']);
tline= fgetl(fid);
while ischar(tline)
    try
    if strcmp(tline(1:20),'TimeStamp of the end')
        tline=tline(59:end);
        EndTimestamp=sscanf(tline,'%d',1);
    end
    catch end
    try
    if strcmp(tline(1:22),'TimeStamp of the start')
        tline=tline(61:end);
        StartTimestamp=sscanf(tline,'%d',1);
    end
    catch
    end
    try
    if strcmp(tline(1:9),'Number of')
        tline=tline(31:end);
        ChanNum=sscanf(tline,'%d',1);
    end
    catch end
    try
    if strcmp(tline(1:9),'File size')
        tline=tline(21:end);
        DatSize=sscanf(tline,'%lu',1);
    end
    catch end

    tline= fgetl(fid);
end
fclose(fid);

%Calculate file length from dat file sample rat
DatLength=DatSize/(ChanNum*2*20) %Dat file size in ms
TspLength=EndTimestamp-StartTimestamp
keyboard
%interpolate to 1 kHz
tspdelta = tspdata(:,1) - tspdata(1,1);
tspdata(:,1) = tspdata(1,1) + tspdelta * DatLength/TspLength;
%interpTsp1kHz(:,1)=tspdata(1,1):tspdata(end,1);
interpTsp1kHz(:,1)=tspdata(1,1):round(tspdata(end,1));
findgood1 = find(tspdata(:,2)>0); if length(findgood1)<2; findgood1 = find(tspdata(:,2)); end   % Color 1
findgood2 = find(tspdata(:,4)>0); if length(findgood2)<2; findgood2 = find(tspdata(:,4)); end   % Color 2
findgood3 = find(tspdata(:,6)>0); if length(findgood3)<2; findgood3 = find(tspdata(:,6)); end   % Color 3

interpTsp1kHz(:,2:3)=floor(interp1(tspdata(findgood1,1),tspdata(findgood1,2:3),interpTsp1kHz(:,1))); % Color1-XY
interpTsp1kHz(:,4:5)=floor(interp1(tspdata(findgood2,1),tspdata(findgood2,4:5),interpTsp1kHz(:,1))); % Color2-XY
interpTsp1kHz(:,6:7)=floor(interp1(tspdata(findgood3,1),tspdata(findgood3,6:7),interpTsp1kHz(:,1))); % Color3-XY


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %align the beginning
    if (StartTimestamp<tspdata(1,1)) 
        i=1:tspdata(1,1)-StartTimestamp;
        tspnew(i,1)=StartTimestamp+i-1;
        tspnew(i,2:7)=-1;
        tspnew=[tspnew; interpTsp1kHz];
    else
        tspnew=interpTsp1kHz(StartTimestamp-tspdata(1,1)+1:end,:);
    end
    clear i

    %adjust the end
    if (TspLength<size(tspnew,1)) 
        tspnew=tspnew(1:TspLength,:);
    else
        %for i=1:TspLength-size(tspnew,1);
        %tspnew=[tspnew; [tspnew(end,1)+1 repmat(-1,1,6)]];
        %end
        addlength = TspLength-size(tspnew,1);
        tspnew = [tspnew; [[(tspnew(end,1)+1):1:(tspnew(end,1)+addlength)]' repmat(-1,addlength,6)]];
    end
    clear i

    clear interpTsp1kHz

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interpTsp(:,1)=tspnew(1,1):size(tspnew,1)/DatLength:tspnew(end,1);
find1 = find(tspnew(:,2)>0); if length(find1)<2; find1 = find(tspnew(:,2)); end   % Color 1
find2 = find(tspnew(:,4)>0); if length(find2)<2; find2 = find(tspnew(:,4)); end   % Color 2
find3 = find(tspnew(:,6)>0); if length(find3)<2; find3 = find(tspnew(:,6)); end   % Color 3

interpTsp(:,2:3) = floor(interp1(tspnew(find1,1),tspnew(find1,2:3),interpTsp(:,1)));
interpTsp(:,4:5) = floor(interp1(tspnew(find2,1),tspnew(find2,4:5),interpTsp(:,1)));
interpTsp(:,6:7) = floor(interp1(tspnew(find3,1),tspnew(find3,6:7),interpTsp(:,1)));
interpTsp(:,8)   = (-1) * ones(size(interpTsp(:,1))); % all -1

interpTsp(isnan(interpTsp))=-1; % extrapolation ponts = NaN --> -1

%fid=fopen([file '.whl1k'],'w');
%for i=1:size(interpTsp,1)
%    fprintf (fid,'%i %i\n',interpTsp(i,2:3));
%end
%fclose(fid);

% Here we don't save the whl1k file - useless ?

%whl1k = interpTsp(:,mycolor);
%save([file '.whl1k'],'whl1k','-ascii')

disp(['Samples in .dat file per channel: ' int2str(DatLength)]);
disp(['Lines in .whl1k file:' int2str(size(interpTsp,1))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%resample to 39.0625 Hz

%whldata(:,1)=floor(resample(interpTsp(:,2),390625,10000000));
%whldata(:,2)=floor(resample(interpTsp(:,3),390625,10000000));
whltime = [interpTsp(1,1):(1000/39.0625):interpTsp(end,1)]';
find1 = find(interpTsp(:,2)>0); if length(find1)<2; find1 = find(interpTsp(:,2)); end   % Color 1
find2 = find(interpTsp(:,4)>0); if length(find2)<2; find2 = find(interpTsp(:,4)); end   % Color 2
find3 = find(interpTsp(:,6)>0); if length(find3)<2; find3 = find(interpTsp(:,6)); end   % Color 3

whldata(:,1)   = whltime;
whldata(:,2:3) = floor(interp1(interpTsp(find1,1),interpTsp(find1,2:3),whltime));
whldata(:,4:5) = floor(interp1(interpTsp(find2,1),interpTsp(find2,4:5),whltime));
whldata(:,6:7) = floor(interp1(interpTsp(find3,1),interpTsp(find3,6:7),whltime));
whldata(:,8)   = (-1) * ones(size(whldata(:,1))); % all -1

whldata(isnan(whldata))=-1; % extrapolation ponts = NaN --> -1

%save it
%fid=fopen([file '.whl'],'w');
%for i=1:size(whldata,1)
%    fprintf (fid,'%i %i\n',whldata(i,:));
%end
%fclose(fid);

whl = whldata(:,mycolor);
save([file '.whl'],'whl','-ascii')



end
