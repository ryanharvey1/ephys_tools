ses = 1; %only look at session 1

for i = 1:length(groupid_place)
    
    temp=load(['F:\ClarkP30_Recordings\ProcessedData\',groupid_place{i,1}],...
        'events','frames','spikesID','Spikes','samplerate','ratemap','maze_size_cm');
    
    cell=find(contains(temp.spikesID.TetrodeNum,groupid_place{i,2}) & ismember(temp.spikesID.CellNum,str2double(groupid_place{i,3})))';
    
    if size(temp.events,2) == 4
        RateMap1 = temp.ratemap{cell,ses};
        RateMap2 = temp.ratemap{cell,ses+1};
        RateMap1(isnan(RateMap1))=0; RateMap1(isinf(RateMap1))=0;
        RateMap1 = padarray(RateMap1,[3 3],0,'both');
        RateMap2(isnan(RateMap2))=0; RateMap2(isinf(RateMap2))=0;
        RateMap2 = padarray(RateMap2,[3 3],0,'both');
        similarity(i,1) = corr2(RateMap1,RateMap2);
    else
        similarity(i,1) = nan; 
    end
    
    clear RateMap1 RateMap2 temp cell
%     fig=figure('Name',id{i,1},'NumberTitle','off');
%     postprocessFigures.ratemaps_2d(fig,temp.ratemap{cell,ses},group_place(i,:,ses),varnames)
%    
%     export_fig('d:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Projects\ClarkP30_Ephys\Analysis\Place_data\Figures','-eps'
end

% Phase precession (spike on phase) plots for Tg+ and WT rats 

