%temporal analysis of OF paths LB July 2019
%Purpose of this script is to evaluate the paths of animals on an open
%field over time. 

%Get data table

cd('D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data')
load('params_V17.mat') %%loads table created from OF_preprocess.
clearvars -except params
param_idx=params.subID; %serves as index for getParam
fr=30;
bin=120; %desired time bin in seconds
 params.distances{1}=[];
for j=1:size(param_idx,1) %loop within subjects
    back=params.backCM{j};
    ts=params.ts{j};
    time_interval=1:bin*fr:size(back,1); 
    time_interval=time_interval'; %creates column vector of time interval starts for indexing below
    
    for i=1:size(time_interval,1) %loop within time intervals. 
        if i==size(time_interval,1)
            temp=back(time_interval(i,1):size(back,1),:); %set last time interval to last time bin to end of session
        else
            temp=back(time_interval(i,1):time_interval(i+1)-1,:);
        end
        for hb=1:size(params.HBcenter{j},2)
            distances(i,hb) = nanmean(sqrt(sum(bsxfun(@minus, temp, [params.HBcenter{j}{1,hb}(:,1),params.HBcenter{j}{1,hb}(:,2)]).^2,2)));
        end
    end
    params.distances{j}=distances; clear distances
end
