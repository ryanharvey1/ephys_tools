function ok=cfsession(action,fn)
%CFSESSION Clear, load, or save a Curve Fitting session

%   $Revision: 1.15.2.12 $  $Date: 2008/05/01 20:12:37 $
%   Copyright 2001-2008 The MathWorks, Inc.

if nargin<2
   fn = '';
end

% Turn off legend before loading; loaded session may turn it back on
oldleg = cfgetset('showlegend');
if isequal(oldleg,'on') && isequal(action,'load')
   cfgetset('showlegend','off');
end

% Wrap real function with a try/catch block so we can catch/display errors
ok = false;
try
   ok = cfsession1(action,fn,oldleg);
catch e
   uiwait(errordlg(sprintf('Unexpected error occurred: %s', e.message ),...
                   'Curve Fitting Tool','modal'))
end

% Session is no longer dirty (in need of saving) if we succeeded
if ok
   cfgetset('dirty',false);
end



% --------- Helper function does the real work
function ok=cfsession1(action,fn,oldleg)

import com.mathworks.toolbox.curvefit.*;

ok = 1;
ftypenow = 'Curve Fitting session';
versionnow = 2;       % the most current version
allversions = [1 2];  % a list of all supported versions

% version 1 was the original
%         2 added the 'guistate' field, but should be compatible with 1

% Variables we save, and the number required to be in a saved file
varnames = {'ftype' 'version' 'allds' 'dsinfo' 'allfits' 'allcfits' ...
            'fitinfo' 'customlibrary' 'outliers' 'guistate'};
nrequired = 8;

switch(action)
 % ---- Save current Curve Fitting tool session
 case 'save'
   % Get file name to use, remember the directory name
   olddir = cfgetset('dirname');
   filespec = [olddir '*.cfit'];
   if nargin<2 || isempty(fn)
      [fn,pn] = uiputfile(filespec,'Save Session');
      if isequal(fn,0) || isequal(pn,0)
         ok = 0;
         return
      end
      if ~ismember('.',fn)
         fn = [fn '.cfit'];
      end
      cfgetset('dirname',pn);
      fn = [pn fn];
   end

   ftype = ftypenow;           %#ok variable saved in varnames 
   version = versionnow;       %#ok variable saved in varnames
  
   % Get all M data set object instances and some properties
   allds = cfgetalldatasets;
   nds = length(allds);
   dsinfo = cell(nds,1);       %#ok variable saved in varnames
   for j=1:length(allds)
      % Save all datasets and their plot info
      dj = allds{j};
      dsinfo{j} = dj.plot;      
   end
   
   % Get all M fit object instances and some properties
   allfits = cfgetallfits;
   nfits = length(allfits);
   allcfits = cell(size(allfits));  %#ok variable saved in varnames
   fitinfo = cell(size(allfits));   %#ok variable saved in varnames
   for j=1:nfits
      % Save cfit objects separately, remember line properties
      fj = allfits{j};
      fitinfo{j} = fj.plot;
      allcfits{j} = fj.fit;
   end
   
   % get the user-defined custom equations
   customlibrary=cfgetset('customLibrary'); %#ok variable saved in varnames
   
   % Get the outliers (excluded sets)
   outdb = getoutlierdb;
   outliers = find(outdb);
   outliers(outliers==outdb) = [];

   % Get information about the current gui state
   guistate.residptype = cfgetset('residptype');
   guistate.showlegend = cfgetset('showlegend');
   if isequal(guistate.showlegend,'on')
      cffig = cfgetset('cffig');
      ax = findall(cffig,'Type','axes','Tag','main');
      if isscalar(ax)
         legh = legend(ax);
         if isscalar(legh)
            guistate.mainrelpos = getrelativelegendposition(cffig,ax,legh);
            guistate.leginfo = cfgetlegendinfo(legh);
         end
      end
      ax = findall(cffig,'Type','axes','Tag','resid');
      if isscalar(ax)
         legh = legend(ax);
         if isscalar(legh)
            guistate.residrelpos = getrelativelegendposition(cffig,ax,legh);
            guistate.rleginfo = cfgetlegendinfo(legh);
         end
      end
   end

   % Select a file and save the session variables
   try
      save(fn, varnames{:}, '-mat');
   catch e
      uiwait(errordlg(sprintf('Error saving session file: %s', e.message ),...
                      'Save Error','modal'))
      ok = 0;
      return
   end
   ok = 1;
   
 % ---- Load new session into Curve Fitting tool
 case 'load'
   % Get file name and load from it, remember the directory name
   olddir = cfgetset('dirname');
   filespec = [olddir '*.cfit'];
   
   if nargin<2 || isempty(fn)
      [fn,pn] = uigetfile(filespec,'Load Session');
      if isequal(fn,0) || isequal(pn,0)
         return
      end
      if ~ismember('.',fn)
         fn = [fn '.cfit'];
      end
      cfgetset('dirname',pn);
      fn = [pn fn];
   end
   
   % Clear current session
   cfsession('clear');
   com.mathworks.toolbox.curvefit.DataSetsManager.getDataSetsManager.turnOffUDDListener;
   try
      s = load('-mat',fn);
   catch e
      uiwait(errordlg(sprintf('Error loading session file: %s', e.message ),...
                      'Load Error','modal'))
   	  com.mathworks.toolbox.curvefit.DataSetsManager.getDataSetsManager.turnOnUDDListener;
      return
   end
   com.mathworks.toolbox.curvefit.DataSetsManager.getDataSetsManager.turnOnUDDListener;
   
   for j=1:nrequired
      if ~isfield(s,varnames{j})
         uiwait(errordlg('Not a valid Curve Fitting session file',...
                      'File Invalid','modal'))
         return
      end
   end
   if ~isequal(s.ftype,ftypenow)
      uiwait(errordlg('Not a valid Curve Fitting session file',...
                      'File Invalid','modal'))
      return
   end

   if ~ismember(s.version,allversions)
      uiwait(errordlg('Bad version number in Curve Fitting session file',...
                      'Invalid Version','modal'))
      return
   end

   % Plot datasets that are flagged for plotting
   for j=1:length(s.allds);
      dj = s.allds{j};
      dj.line = [];
      updatelim(dj);
      if s.dsinfo{j}
         if dj.plot
            % if plot flag already set, need to trigger listener directly
            cfmplot(dj);
         else
            % otherwise setting the flag will trigger the listener
            dj.plot = 1;
         end
      else
         dj.plot = 0;
      end
      issmooth = java.lang.Boolean(~isempty(dj.Source));
      dslength = java.lang.Double(dj.xlength);
      com.mathworks.toolbox.curvefit.DataSetsManager.getDataSetsManager.addDataSet(...
            java(dj),dj.name,dslength,issmooth);
   end
   
   % Load the outliers (excluded sets)
   if ~isempty(s.outliers)
      outdb = getoutlierdb;
      for j=1:numel(s.outliers)
         % add it to the list of outliers
         connect(s.outliers(j),outdb,'up');
      end
      com.mathworks.toolbox.curvefit.OutliersManager.getOutliersManager.init;
   end

   % Fix up fit objects
   for j=1:length(s.allfits)
      fj = s.allfits{j};

      % Restore all cfit object instances
      fj.fit = s.allcfits{j};
   
      % Restore all dataset handles
      dsname = fj.dataset;
      for k=1:length(s.allds)
         if isequal(dsname,s.allds{k}.name)
            fj.dshandle = s.allds{k};
            break;
         end
      end

      % Plot if so flagged
      if s.fitinfo{j}
         fj.line = [];
         fj.rline = [];
         fj.plot = 1;
      end
      isgood = java.lang.Boolean(fj.isgood);
      com.mathworks.toolbox.curvefit.FitsManager.getFitsManager.addFit(...
            java(fj),fj.name,isgood,fj.outlier,dsname);
   end
   
   % Restore the user-defined custom equations
   cfgetset('customlibrary',s.customlibrary);
   if ~isempty(s.customlibrary)
      names=s.customlibrary.names;
      customEquationList=CustomEquationList.getInstance;
      for i=1:length(names)
         customEquationList.addEquation(names{i});
      end
   end

   % Restore the gui state.  Create default state for old saved sessions.
   if ~isfield(s,'guistate')
      s.guistate.residptype = cfgetset('residptype');
      s.guistate.showlegend = oldleg;
   end

   % Residuals
   if isfield(s.guistate,'residptype')
      cftool('toggleresidplot',s.guistate.residptype);
   end
   
   % Legends and their positions
   if ~isfield(s.guistate,'showlegend')
      s.guistate.showlegend = oldleg;
   end
   if ~isfield(s.guistate,'leginfo')
       s.guistate.leginfo = {};
   end
   if ~isfield(s.guistate,'rleginfo')
       s.guistate.rleginfo = {};
   end
   cftool('togglelegend',s.guistate.showlegend,s.guistate.leginfo,...
                         s.guistate.rleginfo);
   if isequal(s.guistate.showlegend,'on')
      cffig = cfgetset('cffig');
      if isfield(s.guistate,'mainrelpos')
         ax = findall(cffig,'Type','axes','Tag','main');
         if isscalar(ax)
            setrelativelegendposition(s.guistate.mainrelpos,cffig,ax);
         end
      end
      if isfield(s.guistate,'residrelpos')
         ax = findall(cffig,'Type','axes','Tag','resid');
         if isscalar(ax)
            setrelativelegendposition(s.guistate.residrelpos,cffig,ax);
         end
      end
   end
   
 % ---- Clear current Curve Fitting tool session
 case 'clear'
   % Trigger java listeners to clear all saved java content
   CFToolClearManager.getCFToolClearManager.listenerTrigger;

   % Delete all udd fit object instances
   allfits = cfgetallfits;
   for j=1:length(allfits)
      fj = allfits{j};
      delete(fj);
   end
   
   % Delete all udd data set object instances
   allds = cfgetalldatasets;
   for j=1:length(allds)
      dj = allds{j};
      if dj.plot, dj.plot = 0; end
      delete(dj);
   end

   % Delete all udd outlier object instances
   outdb = getoutlierdb;
   outliers = find(outdb);
   outliers(outliers==outdb) = [];
   delete(outliers);

   managecustom('clear'); %delete all custom equations
   
   %init (behind the scenes) analysis and plot
   awtinvoke('com.mathworks.toolbox.curvefit.Analysis', 'initAnalysis');
   awtinvoke('com.mathworks.toolbox.curvefit.Plotting', 'initPlotting');
end

% ---- Remember x and y limits, useful for selecting plot limits
function updatelim(h)

if isempty(h.x)
   h.xlim = [];
else
   h.xlim = [min(h.x) max(h.x)];
end
if isempty(h.y)
   h.ylim = [];
else
   h.ylim = [min(h.y) max(h.y)];
end

