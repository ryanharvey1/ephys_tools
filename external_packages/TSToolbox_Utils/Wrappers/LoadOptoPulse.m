function [t_st t_en l singleIx doubleIx] = LoadOptoPulse(fbasename)

% USAGE 
%     [t_st t_en l singleIx doubleIx] = LoadOptoPulse(fbasename)
% 
% read fbasename_optoStim.mat and retrive the pulse info (not oscillatory)
%     
% OUTPUT:
%     - t_st / t_en: ts of start/end illumination times
%     - l: illumination duration
%     - singleIx = indices of single stim
%     - doubleIx = indices of double stim
%     
% Adrien Peyrache, 2012

if exist([fbasename '_optoStim.mat'],'file')
    load([fbasename '_optoStim.mat'],'pulseStim');
else
    error('No file')
end

t_st = pulseStim(:,1)*10000/1250;
t_en = (pulseStim(:,1)+pulseStim(:,2))*10000/1250;
l = t_en-t_st;

d = t_en(2:end)-t_st(1:end-1);

dix = find(d<5000);
if any(diff(dix)==1)
    warning('Found some triplets...')
end
doubleIx = [dix dix+1];

singleIx = find(~ismember([1:length(t_st)]',unique([dix;dix+1])));

t_st = ts(t_st);
t_en = ts(t_en);
