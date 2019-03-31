
function Process_ConvertBasler2Pos(fbasename,varargin)

syncDatFile = [fbasename '_digitalin.dat'];
syncSampFq  = 20000;
syncChan    = 1;
syncNbCh    = 1;
posSampFq   = 60; %in Hz

% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help Process_ConvertBasler2Pos'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'syncdatfile',
      syncDatFile = varargin{i+1};
      if ~isa(duration,'char')
        error('Incorrect value for property ''syncDatFile'' (type ''help Process_ConvertBasler2Pos'' for details).');
      end
    case 'syncsampfq',
      syncSampFq = varargin{i+1};
      if ~isa(frequency,'numeric')
        error('Incorrect value for property ''syncsampfq'' (type ''help Process_ConvertBasler2Pos'' for details).');
      end
    case 'syncchan',
      syncChan = varargin{i+1};
      if ~isa(start,'numeric')
        error('Incorrect value for property ''syncChan'' (type ''help Process_ConvertBasler2Pos'' for details).');
      end
		if start < 0, start = 0; end
    case 'syncnbch',
      syncNbCh = varargin{i+1};
      if ~isa(syncNbCh,'numeric')
        error('Incorrect value for property ''syncNbCh'' (type ''help Process_ConvertBasler2Pos'' for details).');
      end
    case 'possampfq',
      posSampFq = varargin{i+1};
      if ~isa(posSampFq,'numeric')
        error('Incorrect value for property ''posSampFq'' (type ''help Process_ConvertBasler2Pos'' for details).');
      end
    
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help Process_ConvertBasler2Pos'' for details).']);
  end
end

% try
    
    dat = LoadBinary(syncDatFile,'nchannels',syncNbCh,'channels',syncChan);
    t = (0:length(dat)-1)'/syncSampFq;
    
    if ~exist([fbasename '.led'])    
        pos = Process_DetectLED(fbasename);
        if isempty(pos)
            warning('Seems like video file is absent or corrupted. Pos file will all NaN at 60Hz');
            timestamps = [0:1/60:(length(dat))/syncSampFq-1/60];
            pos = -ones(length(timestamps),4);
            dlmwrite([fbasename '.pos'],[timestamps(:) pos],'delimiter','\t','precision',9);
            return
        end
    else
        pos = load([fbasename '.led']);
    end
    
    %It's a pull-off, trigger is toward 0
    dat = double(dat<1);
    dPos = find(diff(dat)==1);
    dNeg = find(diff(dat)==-1);
    
    if length(dPos) == length(dNeg)+1
        dPos = dPos(1:end-1);
    elseif length(dPos) > length(dNeg)+1 || length(dPos) < length(dNeg)
        keyboard
    end
    
    %Frame timing is the middle of shuter opening
    frameT = (t(dPos)+t(dNeg))/2;
    
    %The camera keeps on recording a few frames after software stopped
    %recording. So we skip the last frames of the TTL
    if length(frameT)<size(pos,1);
        warning('Too many video frames!')
        pos = pos(1:length(frameT),:);
    elseif length(frameT)>size(pos,1);
        frameT = frameT(1:size(pos,1));
    end
    
    %We now interpolate the data at 60 Hz (or sampling fqcy specified in
    %arguments)
    recDuration = length(dat)/syncSampFq;
    
    timestamps = (0:1/posSampFq:recDuration-1/posSampFq)';
    pos(pos==-1) = NaN;
    newPos = interp1(frameT,pos,timestamps);
    
    dlmwrite([fbasename '.pos'],[timestamps newPos],'delimiter','\t','precision',9);
    
% catch
%    keyboard
% end
end

