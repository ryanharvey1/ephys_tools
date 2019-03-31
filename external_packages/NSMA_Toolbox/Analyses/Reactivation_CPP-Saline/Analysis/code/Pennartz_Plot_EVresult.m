%%%
%%% script to plot result of EV run stored in EV_results cell array
%%% expects a NRats*NSessions cell array named EV_results
%%%

% convert info in cell arrays to numeric arrays
clear ev*

%%load EV_3rats_AllSessions100_2.2Hz_Results

[NRats,NSess] = size(EV_results);

ev_all = [];
ev_rev_all = [];
count = 0;
for irat = 1:NRats
    ev{irat} = [];
    ev_rev{irat} = [];
    for isess = 1:NSess
        if ~isempty(EV_results{irat,isess})
            ev{irat} = [ev{irat}, EV_results{irat,isess}{2}]; 
            ev_rev{irat} = [ev_rev{irat}, EV_results{irat,isess}{3}]; 
            count = count + 1;
            ev_all(count) = ev{irat}(end);
            ev_rev_all(count) = ev_rev{irat}(end);
        end
    end
    ev_avg(irat) = mean(ev{irat});
    ev_std(irat) = std(ev{irat})/sqrt(length(ev{irat}));
    ev_rev_avg(irat) = mean(ev_rev{irat});
    ev_rev_std(irat) = std(ev_rev{irat})/sqrt(length(ev_rev{irat}));
    [h,p,ci,stat] = ttest2(ev{irat},ev_rev{irat})
    tt2{irat} = {h,p,ci,stat};
end
ev_all_avg = mean(ev_all);
ev_all_std = std(ev_all)/sqrt(length(ev_all));
ev_rev_all_avg = mean(ev_rev_all);
ev_rev_all_std = std(ev_rev_all)/sqrt(length(ev_rev_all));
[h,p,ci,stat] = ttest2(ev_all,ev_rev_all)
tt2{NRats+1} = {h,p,ci,stat};

%plot result
figure
barerror([1;2;3;4],...
    [[ev_avg,ev_all_avg]',[ev_rev_avg,ev_rev_all_avg]'],...
    [[ev_std,ev_all_std]',[ev_rev_std,ev_rev_all_std]'],...
    {'Rat 1','Rat 2','Rat 3','All Rats'},...
    {'EV','reverse EV'})



