% Mat = Accumulate(Index, Value, Size);
%
% This function does something MATLAB really should have built in.
%
% it performs a matrix accumulation.  i.e. for each Index, it adds
% the corresponding Value to that element of a matrix.
%
% this is using sparse and converting back to a full matrix.
%
% inputs: Index: n by d array giving the d-dimensional indices
% of the n array elements to accumulate into.
%
% Value: n by 1 array giving the value to add to that element.
% Can be a scalar, in which case that value is added to each element.
% default 1.
%
% Size: d by 1 array giving size of output matrix.  default: max(Index).


function Mat = Accumulate(Index, Value, Size);

if nargin<2
	Value = 1;
end

if nargin<3
	Size = max(Index);
end

n = size(Index,1);
d = size(Index,2);

%if prod(size(Index))==max(size(Index))
if size(Index,2)==1
    % need to make 1d 2d because MATLAB is annoying
    Index = [Index(:), ones(length(Index),1)];
    Size = [Size 1];
    d=2;
end

if length(Value)==1
	Value = repmat(Value,n,1);
end

if length(Value)~=length(Value(:))
    error('Value should be 1d');
end
Value = Value(:);

if length(Value)~=n | length(Value(:))~=n
	error('Value should be same length as Index, or a scalar');
end

Ip = Index';
if any(Index(:)<0) | any(Ip(:)>repmat(Size(:),n,1)) | ~isequal(Index, floor(Index))
	error('Index should be integers no bigger than Size');
end

if length(Size)~=d
	error('Size doesn''t match dimension of Index');
end


% make linear index from subscripts
for i=1:d
	Sub{i} = Index(:,i);
end
LinInd = sub2ind(Size, Sub{:});

Vec = sparse(LinInd, 1, Value, prod(Size), 1);
Mat = reshape(full(Vec), Size);
