function SR = RestrictTheSMatrix(S,t0,t1)
% SR = RestrictTheSMatrix(S,t0,t1)
%
% Restricts the Spike Matrix S to the interval [t0, t1] (in timestamp 0.1 msec units)
% If t0 >= t1, an error is produced.
% 
% PL Feb 2003

% check for constent input arguments
if t0 >= t1
    error('input argument t0 is requred to be smaller than t1!');
end    
    

for i = 1:length(S)
   SR{i} = Restrict(S{i},t0,t1);
end