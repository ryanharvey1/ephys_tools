function tf = cfIsInManualLimitMode
%CFISINMANUALLIMITMODE True if CFTOOL is in "manual limit" mode
%
%   TF = CFISINMANUALLIMITMODE
%
%   This method returns true if CFTOOL is in "manual limit mode". If
%   CFTOOL is in control of these limits, then this method will return
%   false. 
%
%   By "manual limit mode" we mean that the user is in explicit control
%   of the axes limits. The user has control of the axes limits if
% 
%    either   * zoom is on
%    or       * pan is on
%    or       * the axis limit control is on 
%
%   Note that the axes may have the X- and Y-Limit Modes set to 'auto'
%   or 'manual', but CFISINAUTOMATICLIMITMODE reflects whether CFTOOL
%   has control of limits or if it has ceded control to the user.

%   Copyright 2007-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $    $Date: 2008/05/01 20:12:33 $ 

cffig = cfgetset( 'cffig' );

% Get hold of the various modes that the user might have enabled
hZoom = zoom( cffig );
hPan  = pan(  cffig );

% We need to check zoom, pan or the axis limit control are on.
isOn = @( h ) isequal( 'on', h );
tf = ...
    isOn( hZoom.Enable ) || ...
    isOn(  hPan.Enable ) || ...
    isOn( cfgetset( 'showaxlimctrl' ) );
