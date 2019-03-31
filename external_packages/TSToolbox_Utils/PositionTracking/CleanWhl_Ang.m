function [cleanAng,epCleaned] = CleanWhl_Ang(ang,ep)

st = Start(ep,'s');
en = End(ep,'s');
daSt = Data(Restrict(ang,st));
daEn = Data(Restrict(ang,en));

dt = median(diff(Range(ang)));

addTime = [];
addAng = [];
epCleaned = ep;

for ii = 1 : length(st) -1;
    
    if st(ii+1) - en(ii)<5
        
        da1 = daEn(ii);
        da2 = daSt(ii+1);
        
        diffAng = da2-da1;
        if abs(diffAng)>pi;
            diffAng = mod(diffAng - 2*pi,2*pi);
        end
        
        if abs(diffAng) < pi/3;
            
            t = [en(ii):dt:st(ii+1)]';
            newAng = interp1([en(ii) st(ii+1)]',[da1 da1+diffAng]',t');
            newAng = mod(newAng,2*pi);
            addAng = [addAng;newAng(:)];
            addTime = [addTime;t(:)];
            epCleaned = union(epCleaned,intervalSet(en(ii),st(ii+1)));
            
            if 0
               figure(1),clf
               plot(Range(ang),Data(ang),'b')
               hold on
               plot(t,newAng,'r')
               
               keyboard
            end
            
        end
        
    end
    
end


if 1
   figure(1),clf
   plot(Range(ang),Data(ang),'b')
   hold on
   plot(addTime,addAng,'r')

end

rg = [Range(ang);addTime];
cleanAng = [Data(ang);addAng];

[rg,ix] = sort(rg,'ascend');
cleanAng = cleanAng(ix);

[rg,ix] = unique(rg);
cleanAng = cleanAng(ix);

cleanAng = tsd(rg,cleanAng);
epCleaned = mergeCloseIntervals(epCleaned,0.001);