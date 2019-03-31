function  FITTYPE_OBJ_A_ = getcoeffmatrix(FITTYPE_OBJ_,varargin)  
% GETCOEFFMATRIX Compute linear coefficients matrix.

%   Copyright 2001-2010 The MathWorks, Inc.
%   $Revision: 1.7.2.6 $  $Date: 2010/10/08 16:37:18 $

if ~isequal( 'custom', category( FITTYPE_OBJ_ ) )
    error(message('curvefit:fittype:getcoeffmatrix:NotCustom'));
end
            
CELL_A_ = linearexprs(FITTYPE_OBJ_);
FITTYPE_OBJ_XDATA_ = varargin{end};
varargin(end) = [];
if ~isempty(varargin)
    FITTYPE_PROBPARAMS = varargin;
else
    FITTYPE_PROBPARAMS = {}; 
end
FITTYPE_OBJ_A_ = zeros(length(FITTYPE_OBJ_XDATA_),length(CELL_A_));
FITTYPE_OBJ_NUMCOEFF = numcoeffs(FITTYPE_OBJ_);
FITTYPE_INPUTS_(FITTYPE_OBJ_NUMCOEFF+1:FITTYPE_OBJ_NUMCOEFF+length(FITTYPE_PROBPARAMS)) = FITTYPE_PROBPARAMS; % syntax from fittype/feval so assignData works
FITTYPE_INPUTS_{FITTYPE_OBJ_NUMCOEFF+length(FITTYPE_PROBPARAMS)+1} = FITTYPE_OBJ_XDATA_; % syntax from fittype/feval so assignData works on prob param

if (isempty(FITTYPE_OBJ_.Aexpr))
    FITTYPE_OBJ_A_ = [];
else
    eval(FITTYPE_OBJ_.assignData);
    eval(FITTYPE_OBJ_.assignProb);

    for FITTYPE_OBJ_I_ = 1:length(CELL_A_)
        try
           FITTYPE_OBJ_A_(:,FITTYPE_OBJ_I_) = eval(CELL_A_{FITTYPE_OBJ_I_});
        catch e
           error(message('curvefit:getcoeffmatrix:linearTermError', CELL_A_{ FITTYPE_OBJ_I_ }, e.message));
        end
    end
end
