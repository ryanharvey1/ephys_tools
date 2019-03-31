%function lHandles  = Lines(x,y,col,style,width)
% plots vertical lines at x coortinated
% from y(1) to y(2) , sets color col
function h = Lines(varargin)
ax = axis;
[x, y,col,style,width] = DefaultArgs(varargin,{[],[], 'r','-',1});
if size(col,1)>size(col,2)
    col = col';
end
if isempty(y)
    if isempty(y); y=ax(3:4); end
    nl= length(x);
    x = x(:)'; y =y(:);
    x = repmat(x,2,1);
    y = repmat(y,1,nl);
    h= line(x,y);
elseif isempty(x)
    if isempty(x) x=ax(1:2); end
    nl= length(y);
    x = x(:); y =y(:)';
    x = repmat(x,1,nl);
    y = repmat(y,2,1);
    h= line(x,y);
else
    if length(x)==2
        nl= length(y);
        x = x(:); y =y(:)';
        x = repmat(x,1,nl);
        y = repmat(y,2,1);
        h= line(x,y);
    else
        nl= length(x);
        x = x(:)'; y =y(:);
        x = repmat(x,2,1);
        y = repmat(y,1,nl);
        h= line(x,y);
    end
end

for i=1:nl
    if isstr(col) 
        mycol = col;
    elseif iscell(col) & length(col)==nl
        mycol = col{i};
    elseif  (isnumeric(col) & size(col,2)==3 & size(col,1)==nl)
        mycol = col(i,:);
    elseif  (isnumeric(col) & size(col,2)==3 & size(col,1)==1)
        mycol = col;
    elseif (isnumeric(col) & size(col,2)==1 & size(col,1)==nl) | (isnumeric(col) & size(col,1)==1 & size(col,2)==nl)
        ColorOrder = get(gca, 'ColorOrder');
        mycol = ColorOrder(col,:);
    else
        mycol='k';
    end
    
    if isstr(style) 
        myst = style;
    elseif iscell(style) &  length(style)==nl
        myst = style{i};
    elseif isnumeric(style) & length(unique(style))<4
        stysel = {'-','--','-.',':'};
        [ust, dummy,ui] = unique(style);
        myst = stysel{ui(i)};
    else
        myst='-';
    end
    
    if isnumeric(width) 
        if length(width)==nl
            mywid = width(i);
        else
            mywid = width;
        end
    else
        mywid =1;
    end
    set(h(i),'Color',mycol);
    set(h(i),'LineStyle',myst);
    set(h(i),'LineWidth',mywid);
end

