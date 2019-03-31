function cfoptstop(state)
%CFOPTSTOP Stop fitting iterations for Curve Fitting GUI.

%   Copyright 2001-2009 The MathWorks, Inc.
%   $Revision: 1.3.2.2 $  $Date: 2009/10/10 20:03:35 $

cfInterrupt( 'set', state );
