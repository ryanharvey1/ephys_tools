function spd = angSpeed(ang,win)

if length(ang)>1
    ang = ang(:);
    da = diff(ang);
    da(da>pi) = da(da>pi)-2*pi;
    da(da<-pi) = da(da<-pi)+2*pi;

    da = abs(gaussFilt(da,win));
    spd = gaussFilt(da,win);
    spd = [spd(1);spd];
else
    spd = zeros(length(ang),1);
end