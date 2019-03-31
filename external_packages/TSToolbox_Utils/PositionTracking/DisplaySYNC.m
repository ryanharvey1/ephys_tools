function [sync1,sync2] = DisplaySYNC4(sync1,sync2,nSamplesPerScreen)

AFdeletedspots = [];
AFaddedspots = [];

% Align first nSamplesPerScreen samples

if nargin < 3,
    nSamplesPerScreen = 100;
end
nSamplesPerScreen = min([nSamplesPerScreen length(sync1) length(sync2)]);
ratio = (sync1(nSamplesPerScreen)-sync1(1))/(sync2(nSamplesPerScreen)-sync2(1));
shift = sync1(1) - sync2(1)*ratio;

% Automatic preselection of potential spurious sync1 points to delete
%     i.e. points i for which t(i)-t(i-1) < 1/2*mean(t)
% and of potential missing sync1 points to add
%     i.e. points i for which t(i)-t(i-1) > 3/2*mean(t)

[syncToDelete,syncToAdd] = GetSuggestions(sync1,1);

% Display nSamplesPerScreen points at a time

figure(3);hold on;
nSync1 = length(sync1);
nSync2 = length(sync2);
n = min([nSync1 nSync2]);
i = 0;
while 1,
    i = i+1;
    if i > ceil(n/nSamplesPerScreen), break; end
    while 1,
        clf;hold on;
        first = (i-1)*nSamplesPerScreen+1;
        last1 = min([i*nSamplesPerScreen nSync1]);
        plot(sync1(first:last1),zeros(last1-first+1,1),...
            '+','markersize',10,'linestyle','none','color',[1 0 0]);
        finddeleted = find(AFdeletedspots>sync1(first) & AFdeletedspots<sync1(last1));
        if ~isempty(finddeleted),
            plot(AFdeletedspots(finddeleted), zeros(length(finddeleted),1),'x','markersize',10,'linestyle','none','color',[0 0 0]);
        end    
        findadded = find(AFaddedspots>sync1(first) & AFaddedspots<sync1(last1));
        if ~isempty(findadded),
            plot(AFaddedspots(findadded), zeros(length(findadded),1),'+','markersize',10,'linestyle','none','color',[0 1 1]);
        end    
        last2 = min([i*nSamplesPerScreen nSync2]);
        plot(sync2(first:last2)*ratio+shift,0.1*ones(last2-first+1,1),...
            '+','markersize',10,'linestyle','none','color',[0 1 0]);
        set(gca,'ylim',[-.025 .125],'xlim',[min([sync1(first) sync2(first)*ratio+shift]) max([sync1(last1) sync2(last2)*ratio+shift])]);
        last = min([i*nSamplesPerScreen nSync1 nSync2]);
        for j = first:last,
            line([sync1(j) sync2(j)*ratio+shift],[0 .1],'color',[0 0 0]);
        end
        
        % Prompt for action
        zoom on;
        action = 'x';
        while ~isempty(action) & action ~= 'd' & action ~= 'D' & action ~= 'a' & action ~= 'A' & action ~= 'F',
            string3 = ' ''F''+ENTER to Fix the video syncronization automatically,';
            if ~isempty(syncToDelete),
                string1 = ' ''D''+ENTER to ask the program to suggest a point to delete,';
            else
                string1 = '';
            end
            if ~isempty(syncToAdd),
                string2 = ' ''A''+ENTER to ask the program to suggest a point to add,';
            else
                string2 = '';
            end
            action = input([' Type' ... 
                    string3 ...
                ' ''d''+ENTER to manually delete a point,' ...
                    string1 ...
                ' ''a''+ENTER to manually add a point,' ...
                    string2 ...  
                ' or just hit ENTER to show the next subset of points...'],'s');
        end
        if isempty(action), break; end
        
        if action == 'd' | (action == 'D' & ~isempty(syncToDelete)),
            
            % Delete point
            
            if action == 'd',
                % Manual selection
                display('   In the figure window, click on the point you wish to delete and hit ENTER');
                [x,y] = ginput;x = x(1);y = y(1);
            end
            if (action == 'd' & y < 0.05) | action == 'D',
                % Correct sync 1
                if action == 'd',
                    % Select SYNC point closest to where the user clicked
                    j = find(sync1 > x);if ~isempty(j), j = j(1); end
                    if j ~= 1 & abs(sync1(j-1)-x) < abs(sync1(j)-x), j = j-1;
                    elseif isempty(j), j = length(sync1); end
                else
                    % Automatic selection of the next potential spurious SYNC point
                    j = syncToDelete(1);
                    zoomIn(sync1(j),0);
                end
                % Plot and ask for confirmation
                p = plot(sync1(j),0,'o','color',[1 0 0],'markersize',16);
                key = '';
                while isempty(key) | (key ~= 'y' & key ~= 'n'),
                    key = input('   Confirm deletion (y/n/z=zoom)? ','s');
                    if ~isempty(key) & key == 'z', zoomIn(sync1(j),0); end
                end
                delete(p);
                if ~isempty(key) & lower(key(1)) == 'y',
                    % Delete
                    sync1(j) = [];
                end
                % Update automatic preselections
                [syncToDelete,syncToAdd] = GetSuggestions(sync1,max([j-1 1])),
            else
                % Correct sync 2 (generally useless because .dat files are reliable)
                j = find(sync2*ratio+shift > x);if ~isempty(j), j = j(1); end
                if j ~= 1 & abs(sync2(j-1)*ratio+shift-x) < abs(sync2(j)*ratio+shift-x), j = j-1;
                elseif isempty(j), j = length(sync2); end
                p = plot(sync2(j)*ratio+shift,0.1,'o','color',[1 0 0],'markersize',16);
                key = input('   Confirm deletion (y/n)?','s');
                delete(p);
                if ~isempty(key) & lower(key(1)) == 'y', sync2(j) = []; end
            end
            
        elseif action == 'a' | action == 'A',
            
            % Add point
            if action == 'a',
                % Manual selection
                display('   In the figure window, click where you wish to add a point and hit ENTER');
                [x,y] = ginput;x = x(1);y = y(1);
            end
            if (action == 'a' & y < 0.05) | (action == 'A' & ~isempty(syncToAdd)),
                % Correct sync 1
                if action == 'a',
                    % Select SYNC point closest to where the user clicked
                    j = find(sync1 > x);if ~isempty(j), j = j(1); end
                    if j == 1,
                        x = sync1(1) - (sync1(2) - sync1(1));
                    elseif isempty(j),
                        j = length(sync1) + 1;
                        x = sync1(end) + (sync1(2) - sync1(1));
                    else
                        x = (sync1(j-1) + median(diff(sync1(1:end))));
                    end
                else
                    % Automatic selection of the next potential missing SYNC point
                    j = syncToAdd(1);
                    x = (sync1(j-1) + median(diff(sync1(1:end))));
                    zoomIn(x,0);
                end
                % Plot and ask for confirmation
                p = plot(x,0,'o','color',[1 0 0],'markersize',16);
                key = '';
                while isempty(key) | (key ~= 'y' & key ~= 'n'),
                    key = input('   Confirm addition (y/n/z=zoom)? ','s');
                    if ~isempty(key) & key == 'z', zoomIn(sync1(j),0); end
                end
                delete(p);
                if ~isempty(key) & lower(key(1)) == 'y',
                    % Add point
                    sync1 = [sync1(1:j-1);x;sync1(j:end)];
                end
                % Update automatic preselections
                [syncToDelete,syncToAdd] = GetSuggestions(sync1,max([j-1 1])),
            else
                % Correct sync 2 (generally useless because .dat files are reliable)
                j = find(sync2*ratio+shift > x);if ~isempty(j), j = j(1); end
                if j == 1,
                    x = (sync2(1) - (sync2(2) - sync2(1)))*ratio+shift;
                elseif isempty(j),
                    j = length(sync2) + 1;
                    x = (sync2(end) + (sync2(2) - sync2(1)))*ratio+shift;
                else
                    x = (sync2(j-1) + sync2(j))/2 *ratio+shift;
                end
                p = plot(x,0,'o','color',[1 0 0],'markersize',16);
                key = input('   Confirm addition (y/n)?','s');
                delete(p);
                if ~isempty(key) & lower(key(1)) == 'y', sync2 = [sync2(1:j-1);(x-shift)/ratio;sync2(j:end)]; end
            end
        elseif action == 'F',
            fprintf('\n\nAUTOFIX iteratively deletes and adds video flashes when adjacent flashes are closer or farther apart than the median. ');
            fprintf('Deleted points are displayed in black. Added points are displayed in cyan.\n');
            
            % add points as necesary beginning with point 2 because point 1 is usually a bit off
            AFdeletedspots = [];
            AFaddedspots = [];
            
            MedSync1Dif = median(diff(sync1(1:end)));
            tolerance = 0.95 ;
            if (median(diff(sync1(2:5))) < median(diff(sync1(1:end)))*tolerance |  median(diff(sync1(2:5))) > median(diff(sync1(1:end)))/tolerance)
                fprintf('\n\n******** The median number of frame between the 2nd-5th spot (%d) is not similar to the median number of frames between all the spots in the file (%d)...',median(diff(sync1(2:5))),median(diff(sync1(1:end))));
                fprintf('i.e there is a problem and you should attempt some manual cleaing before running autofix\n');  
                while 1,
                    goforit = input('do you want to continue? yes/no: ', 's');
                    if strcmp(goforit,'yes') | strcmp(goforit,'no'), break; end
                end
            else
                goforit = 'y';
            end
            if goforit(1) == 'y'              
                [m1,n] = size(sync1);
                % Correct sync 1
                j=2; % skip the first point because it's usually off
                while (j < m1)
                    if (((sync1(j+1)-sync1(j)) >= (MedSync1Dif*tolerance)) & ((sync1(j+1)-sync1(j)) <= (MedSync1Dif/tolerance)))
                        j=j+1;
                    else
                        while ((j < m1) & ((sync1(j+1)-sync1(j)) < (MedSync1Dif*tolerance)))
                            m1=m1-1;              
                            AFdeletedspots = [AFdeletedspots sync1(j+1)];
                            sync1(j+1) = [];           
                        end
                        while ((j < m1) & ((sync1(j+1)-sync1(j)) > (MedSync1Dif/tolerance)))
                            x = sync1(j) + MedSync1Dif;
                            sync1 = [sync1(1:j);x;sync1(j+1:end)];
                            j=j+1;
                            m1=m1+1;
                            m2 = m1;
                            AFaddedspots = [AFaddedspots x];    
                        end
                    end
                end                
                % Update automatic preselections so hitting A or D doesn't confuse the user
                [syncToDelete,syncToAdd] = GetSuggestions(sync1,max([j-1 1])),
            end
            
        elseif strcmp(action,'deleteall'),
            AFdeletedspots = [];
            
            MedSync1Dif = median(diff(sync1(1:end)));
            tolerance = 0.95 ;
            if (median(diff(sync1(2:5))) < median(diff(sync1(1:end)))*tolerance |  median(diff(sync1(2:5))) > median(diff(sync1(1:end)))/tolerance)
                fprintf('\n\n******** The median number of frame between the 2nd-5th spot (%d) is not similar to the median number of frames between all the spots in the file (%d)...',median(diff(sync1(2:5))),median(diff(sync1(1:end))));
                fprintf('i.e there is a problem and you should attempt some manual cleaing before running autofix\n');  
                while 1,
                    goforit = input('do you want to continue? yes/no: ', 's');
                    if strcmp(goforit,'yes') | strcmp(goforit,'no'), break; end
                end
            else
                goforit = 'y';
            end
            if goforit(1) == 'y'
                [m1,n] = size(sync1);
                % Correct sync 1
                j=2; % skip the first point because it's usually off
                while (j < m1)
                    if (((sync1(j+1)-sync1(j)) >= (MedSync1Dif*tolerance)))
                        j=j+1;
                    else
                        while ((j < m1) & ((sync1(j+1)-sync1(j)) < (MedSync1Dif*tolerance)))
                            m1=m1-1;              
                            AFdeletedspots = [AFdeletedspots sync1(j+1)];
                            sync1(j+1) = [];           
                        end
                    end
                end                
                % Update automatic preselections so hitting A or D doesn't confuse the user
                [syncToDelete,syncToAdd] = GetSuggestions(sync1,max([j-1 1])),
            end        
            
         elseif strcmp(action,'addall'),
            AFaddedspots = [];
            
            MedSync1Dif = median(diff(sync1(1:end)));
            tolerance = 0.95 ;
            if (median(diff(sync1(2:5))) < median(diff(sync1(1:end)))*tolerance |  median(diff(sync1(2:5))) > median(diff(sync1(1:end)))/tolerance)
                fprintf('\n\n******** The median number of frame between the 2nd-5th spot (%d) is not similar to the median number of frames between all the spots in the file (%d)...',median(diff(sync1(2:5))),median(diff(sync1(1:end))));
                fprintf('i.e there is a problem and you should attempt some manual cleaing before running autofix\n');  
                while 1,
                    goforit = input('do you want to continue? yes/no: ', 's');
                    if strcmp(goforit,'yes') | strcmp(goforit,'no'), break; end
                end
            else
                goforit = 'y';
            end
            if goforit(1) == 'y'
                [m1,n] = size(sync1);
                % Correct sync 1
                j=2; % skip the first point because it's usually off
                
                while (j < m1)
                    if (((sync1(j+1)-sync1(j)) <= (MedSync1Dif/tolerance)))
                        j=j+1;
                    else
                        while ((j < m1) & ((sync1(j+1)-sync1(j)) > (MedSync1Dif/tolerance)))
                            x = sync1(j) + MedSync1Dif;
                            sync1 = [sync1(1:j);x;sync1(j+1:end)];
                            j=j+1;
                            m1=m1+1;
                            m2 = m1;
                            AFaddedspots = [AFaddedspots x];    
                        end
                    end
                end                
                % Update automatic preselections so hitting A or D doesn't confuse the user
                [syncToDelete,syncToAdd] = GetSuggestions(sync1,max([j-1 1])),
            end        
        end
        
        % Update
        
        nSync1 = length(sync1);
        nSync2 = length(sync2);
        n = min([nSync1 nSync2]);
    end
end

% Helper function: zoom in on coordinates (x,y)

function zoomIn(x,y)

xLim = get(gca,'xLim');
x1 = x - .25*.5*(xLim(2)-xLim(1));
x2 = x + .25*.5*(xLim(2)-xLim(1));
yLim = get(gca,'yLim');
y1 = y - .25*.5*(yLim(2)-yLim(1));
y2 = y + .25*.5*(yLim(2)-yLim(1));
set(gca,'xLim',[x1 x2],'yLim',[y1 y2]);

% Helper function: preselect potential spurious and missing points
% Work only from the index 'start': this is to avoid repeating suggestions
% already discarded by the user

function [syncToDelete,syncToAdd] = GetSuggestions(sync,start)

syncIntervals = diff(sync(start:end));
smallSyncIntervals = find(syncIntervals < 0.95*median(syncIntervals));
syncToDelete = smallSyncIntervals + start;
largeSyncIntervals = find(syncIntervals > median(syncIntervals)/0.95);
syncToAdd = largeSyncIntervals + start;

