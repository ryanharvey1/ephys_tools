function VTfile=FillNaN2( VTfile,column )
%FillNaN2 cycles through angle data, detects sets of NaNs, and replaces them
% with an equally spaced ascending or descending range of numbers from 1
% before to 1 after the set of NaNs.
%   Input: VTfile>> matrix, Column>> column with NaNs
% by Ryan Harvey
for iangle=1:length(VTfile)
    if isnan(VTfile(iangle,column))==1 % IF NAN IS DETECTED
        a=VTfile(iangle-1,column); % ONE ROW ABOVE NAN
        b=VTfile(iangle+1,column); % ONE ROW BELOW NAN
        n=1; % ONE NAN TO BE FILLED
        if isnan(b)==1 % IF ONE ROW BELOW FIRST NAN IS NAN
            for inan=1:length(VTfile) 
                if isnan(VTfile(iangle+inan,column))==1 && isnan(VTfile(iangle+inan+1,column))==0 % COUNTS # OF ROWS OF NANs
                    b=VTfile(iangle+inan+1,column); % # ONE ROW BELOW LAST NAN
                    n=inan+1; % # OF NANS
                    x=linspace(a,b,n); % CREATES RANGE FROM ONE ABOVE TO ONE BELOW NAN
                    for ix=1:length(x) 
                       VTfile(iangle,column)=x(ix); % REPLACE NAN WITH NEW LINESPACE VECTOR
                       iangle=iangle+1; % NEXT ROW
                    end
                    break % RETURN TO MAIN LOOP 
                end
            end 
        elseif isnan(b)==0 
            x=linspace(a,b,n); % CREATES RANGE FROM ONE ABOVE TO ONE BELOW NAN
            for ix=1:length(x)
                VTfile(iangle,column)=x(ix);
                iangle=iangle+1; 
            end
        end
    end 
end 
end

