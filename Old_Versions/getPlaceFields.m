function [fields]=getPlaceFields(varargin)
% USAGE
%
%
% INPUTS
% 
%   ratemap             MxN matrix where M is the number of cells, N is the number
%                       of spatial bins
%   minPeakRate         minimum rate for peak of field [default: 2]
%   minFieldWidth       minimum width of field [default: 2]
%   maxFieldWidth       maximum width of field [default: 30]
%   percentThreshold    percent change between peak rate and start/stop of field
%
%
% OUTPUTS
%
%   fields struct with field data for each cell
%
%
%
% HELP
% This function tries to identify place fields based on a set of heuristics
% 
% written by david tingley, 2017
% Adapted by Ryan Harvey, 2018
debugging_fig=0;

p = inputParser;
addRequired(p,'ratemap',@isnumeric)
addParameter(p,'minPeakRate',2,@isnumeric)
addParameter(p,'minFieldWidth',2,@isnumeric)
addParameter(p,'maxFieldWidth',30,@isnumeric)
addParameter(p,'percentThreshold',.2,@isnumeric)


parse(p,varargin{:})

ratemap = p.Results.ratemap;
minPeakRate = p.Results.minPeakRate;
minFieldWidth = p.Results.minFieldWidth;
maxFieldWidth = p.Results.maxFieldWidth;
threshold = p.Results.percentThreshold; % change between peak rate and start/stop of field


warning off  % findpeaks.m throws warnings if peak isn't found...

for i=1:size(ratemap,1)
    fields{i} = [];
    [pks locs w] = findpeaks(fastrms(ratemap(i,:),5),'minpeakheight',minPeakRate,'MinPeakWidth',minFieldWidth);
    exclude=[];
    for j=1:length(locs)-1
       if min(ratemap(i,locs(j):locs(j+1))) > ((pks(j)+pks(j+1))./2) * threshold
           % exclude fields without a 90 % decrease in rate between peaks
           if pks(j) > pks(j+1)
               exclude = [exclude;j+1];
           elseif pks(j) < pks(j+1)
               exclude = [exclude;j];
           end
       end
    end   
    pks(exclude) = [];
    locs(exclude)=[];
    fieldCount = 1;
    for j=1:length(locs)
        
        Map_Field = ratemap(i,:) > pks(j) * threshold;
        
        start = locs(j);
        stop = locs(j);
        while Map_Field(start) == 1  && start > 1
            start = start-1;
        end
        while Map_Field(stop) == 1  && stop < length(Map_Field) -1
            stop = stop+1;
        end
        if stop - start > minFieldWidth && stop - start < maxFieldWidth
            fields{i}{fieldCount}.start = start;
            fields{i}{fieldCount}.stop = stop;
            fields{i}{fieldCount}.width = stop - start;
            fields{i}{fieldCount}.peakFR = pks(j);
            fields{i}{fieldCount}.peakLoc = locs(j);
            com = start; % calculate center of mass for field
            fields{i}{fieldCount}.COM = fields{i}{fieldCount}.peakLoc;
            while sum(ratemap(i,start:stop)) - sum(ratemap(i,start:com)) > sum(ratemap(i,start:stop))./2
                fields{i}{fieldCount}.COM = com;
                com = com + 1;
            end
            fieldCount = fieldCount + 1;
        end
    end
    
    % if the peak rate is below 1hz...
    if isempty(fields{i})
        [fields{i}{1}.peakFR,fields{i}{1}.peakLoc]=max(ratemap(i,:));
        fields{i}{1}.width=length(ratemap(i,:));
        fields{i}{1}.start=1;
        fields{i}{1}.stop=length(ratemap(i,:));
        fields{i}{1}.COM = fields{i}{1}.peakLoc;
        com=1;
        while sum(ratemap(i,fields{i}{1}.start:fields{i}{1}.stop))...
                - sum(ratemap(i,fields{i}{1}.start:com)) > sum(ratemap(i,fields{i}{1}.start:fields{i}{1}.stop))./2
            fields{i}{1}.COM = com;
            com = com + 1;
        end
    end
    
    % remove fields with the same field boundaries and keep the one with
    % the highest peak rate
    for fie=1:length(fields{i})
        PR(fie)=fields{i}{fie}.peakFR;
        start_(fie)=fields{i}{fie}.start;
        stop_(fie)=fields{i}{fie}.stop;
    end
    fielddets=[start_',stop_',PR'];
    [~,idx]=sort(fielddets(:,3),'descend');
    fielddets=fielddets(idx,:);
    
    [C,ia,~]=unique(fielddets(:,1:2),'rows','stable');
    
    Z=zeros(size(fielddets,1),2);
    Z(ia,:)=C;
    Z=Z(idx,:);
    
    fields_to_delete=find(all(Z==0, 2));
    for f=1:length(fields_to_delete)
        fields{i}{fields_to_delete(f)}=[];
    end
    fields{i}=fields{i}(~cellfun('isempty',fields{i}));

    if debugging_fig
        for f=1:length(fields{i})
            figure;
            plot(ratemap(i,:),'k')
            grid on
            hold on
            plot(fields{i}{f}.start,ratemap(i,fields{i}{f}.start),'*r' )
            text(fields{i}{f}.start,ratemap(fields{i}{f}.start),'start')
            
            plot(fields{i}{f}.stop,ratemap(i,fields{i}{f}.stop),'*r' )
            text(fields{i}{f}.stop,ratemap(i,fields{i}{f}.stop),'stop')
            
            plot(fields{i}{f}.COM,ratemap(i,fields{i}{f}.COM),'*r' )
            text(fields{i}{f}.COM,ratemap(i,fields{i}{f}.COM),'COM')
            
            plot(fields{i}{f}.peakLoc,ratemap(i,fields{i}{f}.peakLoc),'*r' )
            text(fields{i}{f}.peakLoc,ratemap(i,fields{i}{f}.peakLoc),'peak')
        end
    end
end
warning on

