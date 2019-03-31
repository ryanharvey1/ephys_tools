function tbxStruct = demos
%DEMOS   Return HOSA demo information to the MATLAB Demo.

% Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
% $Revision: 1.9 $  $Date: 1999/01/04 21:33:15 $

%     RESTRICTED RIGHTS LEGEND
% Use, duplication, or disclosure by the Government is subject to
% restrictions as set forth in subparagraph (c) (1) (ii) of the 
% Rights in Technical Data and Computer Software clause of DFARS
% 252.227-7013. 
% Manufacturer: United Signals & Systems, Inc., P.O. Box 2374, 
% Culver City, California 90231. 
%
%  This material may be reproduced by or for the U.S. Government pursuant 
%  to the copyright license under the clause at DFARS 252.227-7013. 

if nargout == 0, demo toolbox; return; end

tbxStruct.Name = 'Higher-Order Spectral Analysis';
tbxStruct.Type = 'toolbox';

tbxStruct.Help = { ...
   ' The Higher-Order Spectral Analysis Toolbox is a collection '
   ' of MATLAB functions which implement a variety of advanced '
   ' signal processing algorithms for spectral estimation, '
   ' polyspectral estimation, and computation of time-frequency '
   ' distributions, with applications such as parametric and '
   ' non-parametric blind identification of non-minimum phase '
   ' linear systems, harmonic retrieval in colored noise, time '
   ' delay estimation (TDE) direction of arrival estimation (DOA), ' 
   ' input/output identification of second-order  Volterra models, '
   ' detection of quadratic frequency and phase coupling and '
   ' adaptive AR parameter estimation.  '
   ' '
   ' The toolbox contains mfiles to test for Gaussianity and '
   ' linearity of a time series.  It includes functions for '
   ' estimating cumulants and spectra of orders 2 (spectrum), '
   ' 3 (bispectrum) and 4 (trispectrum).  Bicoherence and cross- '
   ' bicoherence can also be estimated .  It uses cumulants of '
   ' various orders to extract useful information about non- '
   ' Gaussian signals in Gaussian noise '
   ' '
   ' The toolbox is useful for in situations involving the '
   ' so-called "non-problems": linear non-minimum phase systems, '
   ' non-linearities, non-Gaussian signals, non-additive noise, ' 
   ' and non-stationary signals.'
   ' '};

tbxStruct.DemoList = {'Command-line demos', 'hosadem';
                      'Case Studies','ex_case'};
