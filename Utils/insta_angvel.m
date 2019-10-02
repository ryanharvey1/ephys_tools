function anglevel=insta_angvel(theta,samplerate)
% insta_angvel: calc instant angular velocity
%
% Input:
%           theta: angles in degrees
%           samplerate: sample rate
%
% Output:
%           anglevel: angular velocity in angles/sec
%
% Ryan H, Laura B; updated 10/1/2019 
%
normDeg=mod(diff(theta),360);
DiffDeg=[360-normDeg,normDeg];

[~,I]=min(DiffDeg,[],2);

DiffDeg(:,1)=-DiffDeg(:,1);

angvel(I==2,1)=DiffDeg(I==2,2);
angvel(I==1,1)=DiffDeg(I==1,1);

anglevel=angvel*samplerate;

end




