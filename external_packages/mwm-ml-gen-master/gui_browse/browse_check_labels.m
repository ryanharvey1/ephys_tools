function [new_labels,new_indexes] = browse_check_labels(num_tags,labels,indexes)
%BROWSE_CHECK_LABELS modifies a labels cell array according to the number
% of the available tags 

    new_tag = {};
    for i = 1:size(labels,2)
        %get the column
        col = labels(:,i);
        %remove empty cells
        col = col(~cellfun('isempty',col));
        if isempty(col)
            continue;
        end
        %find unique values
        col = unique(col);
        %check if we have any new tag
        for j = 1:length(col)
            if isempty(find(ismember(num_tags,col{j})))
                % keep the new tag
                new_tag = [new_tag,col{j}];
            end
        end
    end
    % if we have new tags ask the user if he wants them removed or conf
    if ~isempty(new_tag);
        str = strcat('New tags were found: ',strjoin(new_tag,','),'. Would you like to delete or configure them?');
        choice = questdlg(str,'New tags','Delete','Configure','Delete');
        if isempty(choice)
            return
        elseif isequal(choice,'Configure');
            msgbox('Press the Configure Tags button','Info','help');
            return
        elseif isequal(choice,'Delete');
            for i = 1:size(labels,2)
                col = labels(:,i);
                for c = 1:length(col)
                    for k = 1:length(new_tag)
                        if isequal(col{c},new_tag{k})
                            col{c} = [];
                        end
                    end
                end
                labels(:,i) = col;
            end
            %remove empty rows and indexes
            emptyCells = cellfun(@isempty,labels);
            rem = [];
            for i = 1:size(emptyCells,1)
                if sum(emptyCells(1,:)) == size(emptyCells,2)
                    rem = [rem,i];
                end
            end
            indexes(rem) = [];
            labels(all(cellfun(@isempty,labels),2),:) = [];
        end
    end
    % if we do not have new tags then labels can be of length num_tags
    if size(labels,2) < length(num_tags)
        labels{1,length(num_tags)} = [];
    elseif size(labels,2) > length(num_tags)
        % in this case just remove columns
        labels = labels(:,1:length(num_tags));
    end

    new_labels = labels;
    new_indexes = indexes;
end

