function cfcshelpviewer(topic, errorname)
% CFCSHELPVIEWER  is a helper file for the Curve Fitting Toolbox 
% CFCSHELPVIEWER Displays help for Curve Fitting TOPIC. If the map file 
% cannot be found, an error is displayed using ERRORNAME

%   Copyright 2001-2008 The MathWorks, Inc.
%   $Revision: 1.4.2.6 $

mapfilename = fullfile( docroot, 'toolbox', 'curvefit', 'csh', 'curvefit_csh.map' );
try
    helpview( mapfilename, topic, 'CSHelpWindow', cfgetset( 'cffig' ) );
    
catch ignore %#ok<NASGU> don't populate the last-error stack
    message = sprintf('Unable to display help for %s\n', errorname);
    errordlg(message);
end
