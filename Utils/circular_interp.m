function spk_phase =  circular_interp(ts,theta,spks)
% interpolates spike times to angle data using sin/cosine 
%input:
%  ts: timestamp vector (length n)
%  theta: head angle in degrees (length n)
%  spks: spike times 
%output: 
%  spk_phase: angles in degrees when spike occured
% get sin and cos then interp spike times
spk_phase = atan2(interp1(ts,sin(deg2rad(theta)),spks,'linear'),interp1(ts,cos(deg2rad(theta)),spks,'linear'));
% wrap to 360
spk_phase = wrapTo360(rad2deg(spk_phase));

end