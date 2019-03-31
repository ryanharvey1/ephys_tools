function CalibrateCamera(varargin)

% CalibrateCamera - Calibrate camera, from pixel to cms
%
%  USAGE
%
%    CalibrateCamera(<options>)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'epoch'       intervalSet objet
%     'diag'        size of the environment diagonal in cms
%     'bon2cms'     bin to cm ratio
%    =========================================================================

% Adrien Peyrache, 2016

fbasename           = extractfbasename(pwd);
[X,Y,GoodRanges]    = LoadPosition_Wrapper(fbasename,0);

diag = 0;
ep = [];
bin2cms = -1;

for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help LoadBinary'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'diag',
      diag = varargin{i+1};
    case 'epoch',
      ep = varargin{i+1};
    case 'bin2cms',
      bin2cms = varargin{i+1};  
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help LoadBinary'' for details).']);
  end
end

if isempty(ep)
    epoch = load('Analysis/BehavEpochs.mat','wakeEp');
    epoch = epoch.wakeEp;
else
    load('Analysis/BehavEpochs',ep)
    eval(['epoch = ' ep ';'])
end

if bin2cms < 0
    epoch   = intersect(epoch,GoodRanges);
    X       = Restrict(X,epoch);
    Y       = Restrict(Y,epoch);

    figure(1),clf
    plot(Data(X),Data(Y))

    done=0;
    while ~done
        fprintf('Draw a line on the figure\n')
        gg = ginput(2);
        lpixel = norm(gg(1,:)-gg(2,:));

        if diag == 0
            diag = input('What is the actual size of this line (in cms)','s');
            diag = str2num(diag);
        end

        bin2cms = diag/lpixel;

        fprintf('Ratio cms/pixel is: %f\n',bin2cms)
        answer = input('Satisfied? [Y/N]','s');
        if strcmp(answer,'Y')
            done=1;
        elseif strcmp(answer,'N')
            fprintf('OK, then let''s restart!\n')
        else
            fprintf('Don''t understand the answer, so we restart\n')
        end
    end
end

if isempty(ep)
    info = {'Bins to cms ratio'};
    varName = '';
else
    info    = {['Bins to cms ratio in ' ep ' epoch']};
    varName = ['_' ep];
end
    
SaveAnalysis(pwd,'GeneralInfo',{bin2cms},{['bin2cms' varName]},info);
