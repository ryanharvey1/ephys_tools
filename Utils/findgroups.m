function [start,ends,ngroups]=findgroups(x)
% findgroups: locates the start and end index of a congigous binary index
% Often used after contigousframes.m to find the locations of groups of contigous frames
% and the corresponding number of groups 
% 
% Input:
%           x: binary index
% Output:
%           start: vector containing start row or col numbers
%           ends: vector containing end row or col numbers
%           ngroups: number of groups within your index
% 
% example:
%       x=[0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 1 1 1 0 0 0 0 1 1 1];
%       [start,ends,ngroups]=findgroups(x)
% 
%       start =
%           8    16    23
% 
%       ends =
%           12    18    25
% 
%       ngroups =
%           3
% 
% Ryan Harvey 2018

% make row vector
x=reshape(x,1,[]);

% pack with zeros
x=[0,x,0];

% locate number of groups
ngroups=sum(diff(x)==1);

% locate start of group
start=find(diff(x)==1);

% locate end of group
x=fliplr(x);
ends=fliplr(abs(find(diff(x)==1)-length(x)))-1;
end

