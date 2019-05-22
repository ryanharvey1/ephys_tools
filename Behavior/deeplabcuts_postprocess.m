% deeplabcuts_postprocess RH May 2019

path_to_files='F:\Users\BClarkLab\Desktop\Maria\WorkDir\TrackedData\lgOF_dlc';

files = dir([path_to_files,'\**\*.csv']);

for i=1:length(files)
    % load header
    fileID = fopen(fullfile(files(i).folder,files(i).name),'r');
    dataArray = textscan(fileID, '%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]',...
        3-2+1, 'Delimiter',',', 'TextType', 'string', 'HeaderLines',...
        2-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    header = [dataArray{1:end-1}];
    clearvars fileID dataArray ans;
    
    % load data
    fileID = fopen(fullfile(files(i).folder,files(i).name),'r');
    dataArray = textscan(fileID, '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]',...
        'Delimiter', ',','TextType', 'string', 'EmptyValue', NaN,...
        'HeaderLines' ,4-1,'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    tsxy= [dataArray{1:end-1}];
    clearvars fileID dataArray ans;
    
    
    % convert ts into seconds
    tsxy(:,1)=0:1/30:tsxy(end,1)/30;
    [tsxy(:,2),tsxy(:,3)]=FixPos(tsxy(:,2),tsxy(:,3),tsxy(:,1));
    [tsxy(:,5),tsxy(:,6)]=FixPos(tsxy(:,5),tsxy(:,6),tsxy(:,1));
    [tsxy(:,8),tsxy(:,9)]=FixPos(tsxy(:,8),tsxy(:,9),tsxy(:,1));
    [tsxy(:,11),tsxy(:,12)]=FixPos(tsxy(:,11),tsxy(:,12),tsxy(:,1));
    
    
    % find single start and end time
    mean_likelihood=mean(tsxy(:,contains(header(2,:),'likelihood')),2);
    ngroups=0;
    while ngroups~=1
        [indexout]=contiguousframes(~[mean_likelihood<mean(mean_likelihood)-2*std(mean_likelihood)],...
            round((100-2).*rand(1,1) + 2));
        [start,ends,ngroups]=findgroups(indexout);
    end
    
    % find head direction
    tsxy(:,14)=wrapTo360(fixNLXangle(rad2deg(atan2(tsxy(:,3)-tsxy(:,6),...
        tsxy(:,2)-tsxy(:,5))),round(0.1667*30)));
    header(1,14)="headdirection";
    header(2,14)="angle";

    % plot
    nose_x = tsxy(start:ends,2);
    nose_y = tsxy(start:ends,3);
    head_x=tsxy(start:ends,5);
    head_y=tsxy(start:ends,6);
    back_x = tsxy(start:ends,8);
    back_y = tsxy(start:ends,9);
    butt_x=tsxy(start:ends,11);
    butt_y=tsxy(start:ends,12);
    
    figure;
    plot(nose_x,nose_y,'.k')
    hold on
    nPlot = plot(NaN,NaN,'ro','MarkerFaceColor','r');
    hPlot = plot(NaN,NaN,'bo','MarkerFaceColor','b');
    bPlot = plot(NaN,NaN,'go','MarkerFaceColor','g');
    buttPlot = plot(NaN,NaN,'co','MarkerFaceColor','c');
    
    
    xlim([min(nose_x) max(nose_x)]);
    ylim([min(nose_y) max(nose_y)]);
    % iterate through each point on line
    for k=1:length(nose_x)
        set(nPlot,'XData',nose_x(k),'YData',nose_y(k));
        set(hPlot,'XData',head_x(k),'YData',head_y(k));
        set(bPlot,'XData',back_x(k),'YData',back_y(k));
        set(buttPlot,'XData',butt_x(k),'YData',butt_y(k));
        
        pause(1/240);
    end
    close
end


