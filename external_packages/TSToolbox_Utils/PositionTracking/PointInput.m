function [x,y,button] = PointInput(varargin)
%function [x,y,button] = PointInput(n,figh, axish)
[n, figh, axish] = DefaultArgs(varargin,{1,gcf,gca});

if n>0
    for i=1:n
        waitforbuttonpress;
        whatbutton = get(figh,'SelectionType');
        mousecoord = get(axish,'CurrentPoint');
        x(i)=mousecoord(1,1);
        y(i) = mousecoord(1,2);

        switch whatbutton
            case 'normal'
                button(i)=1;
            case 'extend'
                button(i)=2;
            case 'alt'
                button(i)=3;
        end
    end
else
    whatbutton = get(figh,'SelectionType');
    mousecoord = get(axish,'CurrentPoint');
    x= mousecoord(1,1);
    y = mousecoord(1,2);

    switch whatbutton
        case 'normal'
            button=1;
        case 'extend'
            button=2;
        case 'alt'
            button=3;
        otherwise 
            fprintf('you pressed button %s\n',whatbutton);
    end
end

