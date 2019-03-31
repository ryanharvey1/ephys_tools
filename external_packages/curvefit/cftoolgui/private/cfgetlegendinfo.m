function leginfo = cfgetlegendinfo(legh)
%CFGETLEGENDINFO Get information about legend so we can re-create it later.

%   $Revision: 1.1.6.1 $  $Date: 2006/06/20 20:04:17 $
%   Copyright 2006 The MathWorks, Inc.

if ~isempty(legh) && ishandle(legh)
    % Record the orientation
    leginfo = {'Orientation', get(legh,'Orientation')};

    % Record other properties that take non-default values
    props = {'EdgeColor'      'DefaultAxesXColor'
             'Color'          'DefaultAxesColor'
             'FontName'       'DefaultTextFontName'
             'FontAngle'      'DefaultTextFontAngle'
             'FontSize'       'DefaultTextFontSize'
             'FontWeight'     'DefaultTextFontWeight'
             'TextColor'      'DefaultTextColor'
             'LineWidth'      'DefaultAxesLineWidth'};
     % "interpreter" is intentionally omitted, as it could affect data set
     % and fit name display
     
     for j=size(props,1):-1:1
         propval = get(legh,props{j,1});
         if isequal(propval, get(0,props{j,2}))
             props(j,:) = [];
         else
             props{j,2} = propval;
         end
     end
     props = props';
     leginfo = [leginfo, props(:)'];
     
     % Include location only if not 'none'
     legloc = get(legh, 'Location');
     if ~isequal(legloc,'none')
         leginfo(end+(1:2)) = {'Location', legloc};
     end
else
    leginfo = {};
end