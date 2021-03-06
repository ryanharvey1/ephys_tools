function data = get_session_idx(data)
% get_session_idx returns session index while compensating for track data 
% which is split into both running directions
%
% Input: Data: ephys_tools data structure
%
% Output: idx: column index for postprocess saving
%
% Ryan H 2020

idx = [];
events=1:length(data.mazetypes);
for event=events
    for i=1:length(data.Spikes)
        % if current event is track and no previous tracks exist
        if contains(data.mazetypes(event),'track','IgnoreCase',true) &&...
                ~any(contains(data.mazetypes(1:event-1),'track','IgnoreCase',true))
            idx = [idx;event;event+1];
            
            % if current event is track and previous tracks exist
        elseif contains(data.mazetypes(event),'track','IgnoreCase',true) &&...
                any(contains(data.mazetypes(1:event-1),'track','IgnoreCase',true))
            n_tracks = sum(contains(data.mazetypes(1:event-1),'track','IgnoreCase',true));
            idx = [idx;event+n_tracks;event+n_tracks+1];
            
            % if current event is not track and no previous tracks exist
        elseif ~contains(data.mazetypes(event),'track','IgnoreCase',true) &&...
                ~any(contains(data.mazetypes(1:event-1),'track','IgnoreCase',true))
            idx = [idx;event];
            
            % if current event is not track and previous tracks exist
        elseif ~contains(data.mazetypes(event),'track','IgnoreCase',true) &&...
                any(contains(data.mazetypes(1:event-1),'track','IgnoreCase',true))
            n_tracks = sum(contains(data.mazetypes(1:event-1),'track','IgnoreCase',true));
            idx = [idx;event+n_tracks];
        end
    end
end
data.session_idx = idx;
end