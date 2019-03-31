function obj = loadobj(obj)
%LOADOBJ Method to post-process FITTYPE objects after loading

%   Copyright 2001-2008 The MathWorks, Inc.
%   $Revision: 1.6.2.2 $  $Date: 2008/10/31 05:56:42 $

% If the loaded version is a structure then we need to convert it to
% an object
if isstruct( obj )
    if obj.version == 1.0
        obj = versionTwoFromVersionOne( obj );
    end
    obj = copyFittypeProperties( fittype, obj );
end

% Restore function handles that had to be removed during save
if ~isequal(category(obj),'custom')
    libname = obj.fType;
    if ~isempty(libname)
        obj = sethandles(libname,obj);
    end
end
end

function obj = versionTwoFromVersionOne( obj )
% VERSIONTWOFROMVERSIONONE -- Convert the old fildnames to the new names
obj = changeFieldname( obj, 'category',        'fCategory' );
obj = changeFieldname( obj, 'constants',       'fConstants' );
obj = changeFieldname( obj, 'feval',           'fFeval' );
obj = changeFieldname( obj, 'fitoptions',      'fFitoptions' );
obj = changeFieldname( obj, 'nonlinearcoeffs', 'fNonlinearcoeffs' );
obj = changeFieldname( obj, 'startpt',         'fStartpt' );
obj = changeFieldname( obj, 'type',            'fType' );
obj = changeFieldname( obj, 'typename',        'fTypename' );
obj.version = 2.0;
end

function s = changeFieldname( s, old, new )
s.(new) = s.(old);
s = rmfield( s, old );
end
