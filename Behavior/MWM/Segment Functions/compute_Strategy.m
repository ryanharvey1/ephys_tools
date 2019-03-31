function [time_quadT,time_quadA,time_quadB,time_quadC] = compute_Strategy(x_smooth,y_smooth,platX,platY,frameRate,r)
%LB 9/2017
%compute_QuadEdge calculates dwell time in quadrants outside an inner
%annulus boundary of a circular maze. 
%this function assumes maze is split into four equal qadrants and coordinates
%are centered around 0,0. Adjacent and opposit quadrants
%are labeled in a clockwise manner from A:C.
%input:
%       -x_smooth: smooth x-coordinates (column vector)
%       -y_smooth: smooth y-coordinates (column vector)
%       -targetQuad: location of target of target quadrant in string form.
%           i.e. 'x,y' indicates target quadrant is within the positive
%           domains of x and y.
%       -frameRate: scalar in frames/s
%       -rOut: radius of pool (cm). 
%       -rIN: radius of inner cannulus (cm)
%       -mkplot: string argument. 'on' will plot coordinates outside of
%       annulus. 'off' does nothing. 
%output:
%       -time_quadT: Total 
%       -time_quadT: Total 
%       -time_quadT: Total 
%       -time_quadT: Total 


%Create function for pool diameter for plots
cosFxnOut = cos(0:2*pi/360:2*pi)*r;
sinFxnOut = sin(0:2*pi/360:2*pi)*r;

%create function for inner annulus
cosFxnIN = cos(0:2*pi/360:2*pi)*r;
sinFxnIN = sin(0:2*pi/360:2*pi)*r;

%create index of points inside inner annulus. 
[in,~]=inpolygon(x_smooth,y_smooth,cosFxnIN,sinFxnIN);

%% computer parameters needed to determine strategies

%time in wall zone 
cosFxnWZ = cos(0:2*pi/length(x_smooth):2*pi)*(r*(.9)); %create circle that is 90% of pool
sinFxnWZ = sin(0:2*pi/length(x_smooth):2*pi)*(r*(.9));
%create index of points inside of annulus. 
[inWZ,~]=inpolygon(x_smooth,y_smooth,cosFxnWZ,sinFxnWZ);
%time in wider wall zone
cosFxnWWZ = cos(0:2*pi/length(x_smooth):2*pi)*(r*(.75)); %create circle that is 75% of pool
sinFxnWWZ = sin(0:2*pi/length(x_smooth):2*pi)*(r*(.75)); %create circle that is 75% of pool
%annuli for chaining max zone 
cosFxnCMax = cos(0:2*pi/length(x_smooth):2*pi)*(r*(.65)); %create circle that is 75% of pool
sinFxnCMax = sin(0:2*pi/length(x_smooth):2*pi)*(r*(.65)); %create circle that is 75% of pool
%annuli for chaining min zone 
cosFxnCMin = cos(0:2*pi/length(x_smooth):2*pi)*(r*(.35)); %create circle that is 75% of pool
sinFxnCMin = sin(0:2*pi/length(x_smooth):2*pi)*(r*(.35)); %create circle that is 75% of pool
%create index outside of annulus. 
[inWWZ,~]=inpolygon(x_smooth,y_smooth,cosFxnWWZ,sinFxnWWZ);
                %mean distance to plat
                %path centroid to plat
                %eccentricity
%time in plat radius (6x platform radius diameter)
r=(rPlat*6); ang=0:0.01:2*pi; xp=r*cos(ang); yp=r*sin(ang);
[inPlatR,~]=inpolygon(x_smooth,y_smooth,(platX+xp),(platY+yp));
                %Detect loop 
                  %Determine Strategy 
                    %Thigmotaxis >35%time in closer wallzone and/or
                    %>65%time in wider wall zone
                    %Chaining >80% time in annulus zone
                    %Focal Search <.35 mean dist. to swim path centroid &
                    %or <.3 mean distand to present goal 
                    %Self orienting: Segment contains a loop 
                    %Perseverance <.35 mean dist. to swim path centroid &
                    %or <.3 mean distand to prior goal (matching to place
                    %only) 
end
