% CCGbias_script
%
% get bias estimate for all 3 epochs
% need sysCCGxxall
%
% KL Hoffman 12/2001

global normal
global RATHC_DATA
s1allzeros = [];
ballzeros = [];
s2allzeros = [];    

disp ('Computing the overall bias for the CCG data of all three epochs')

for isys = 1:length(sysCCGS1all)
    for jsys = 1:length(sysCCGS1all)
        lt = size(sysCCGS1all{isys,jsys}); %length of cellpair vector
        pre = zeros(lt(1),1);
        post = zeros(lt(1),1);
        
        for imat = 1:lt(1)
                 
            pre(imat) = sum(sysCCGS1all{isys,jsys}(imat,31:50));
            post(imat) = sum(sysCCGS1all{isys,jsys}(imat,52:71));
            %if (pre(imat) == 0) | (post(imat) == 0) - only required if you are
            %                                           using some empty tfiles
            %    s1allzeros = [s1allzeros imat];
            %end %if
        end
        
        if (normal)
            S1allbias{isys,jsys} = (post - pre)./ (post + pre);
        else
            S1allbias{isys,jsys} = (post - pre);
        end %if
    end % for jsys
end % for isys
            

for isys = 1:length(sysCCGBall)
    for jsys = 1:length(sysCCGBall)
        lt = size(sysCCGBall{isys,jsys}); %length of cellpair vector
        pre = zeros(lt(1),1);
        post = zeros(lt(1),1);

        for imat = 1: lt(1)
            pre(imat) = sum(sysCCGBall{isys,jsys}(imat,31:50));
            post(imat) = sum(sysCCGBall{isys,jsys}(imat,52:71));
        end
        
        if (normal)
            Ballbias{isys,jsys} = (post - pre)./ (post + pre);
        else
            Ballbias{isys,jsys} = (post - pre);
        end %if
    end % for jsys
end % for isys
    
for isys = 1:length(sysCCGS2all)
    for jsys = 1:length(sysCCGS2all)
        lt = size(sysCCGS2all{isys,jsys}); %length of cellpair vector
        pre = zeros(lt(1),1);
        post = zeros(lt(1),1);
        
        for imat = 1: lt(1)
            pre(imat) = sum(sysCCGS2all{isys,jsys}(imat,31:50));
            post(imat) = sum(sysCCGS2all{isys,jsys}(imat,52:71));
        end
        if (normal)
            S2allbias{isys,jsys} = (post - pre)./ (post + pre);
        else
            S2allbias{isys,jsys} = (post - pre);
        end %if
    end % for jsys
end % for isys


% NOW COMPUTE THE BIAS FOR EACH DATASET 
disp ('Computing the bias for each dataset')
dsetnum = size(sysCCGS1split);
dsetnum = dsetnum(2);
arraynum = length(sysCCGS1split{1});

for isys = 1:arraynum
    for jsys = 1:dsetnum
        lt = size(sysCCGS1split{jsys}{isys}); %length of cellpair vector
        pre = zeros(lt(1),1);
        post = zeros(lt(1),1);
        
        for imat = 1:lt(1)
            pre(imat) = sum(sysCCGS1split{jsys}{isys}(imat,31:50));
            post(imat) = sum(sysCCGS1split{jsys}{isys}(imat,52:71));
            %if (pre(imat) == 0) | (post(imat) == 0)    - only required if there are tfiles
            %                                              without spikes
            %    s1allzeros = [s1allzeros imat];
            %end %if
        end
        
        if (normal)
            S1sessbias{isys,jsys} = (post - pre)./ (post + pre);
        else
            S1sessbias{isys,jsys} = (post - pre);
        end %if
    end % for jsys
end % for isys
            

for isys = 1:arraynum
    for jsys = 1:dsetnum
        lt = size(sysCCGBsplit{jsys}{isys}); %length of cellpair vector
        pre = zeros(lt(1),1);
        post = zeros(lt(1),1);

        for imat = 1: lt(1)
            pre(imat) = sum(sysCCGBsplit{jsys}{isys}(imat,31:50));
            post(imat) = sum(sysCCGBsplit{jsys}{isys}(imat,52:71));
        end
        
        if (normal)
            Bsessbias{isys,jsys} = (post - pre)./ (post + pre);
        else
            Bsessbias{isys,jsys} = (post - pre);
        end %if
    end % for jsys
end % for isys
    
for isys = 1:arraynum
    for jsys = 1:dsetnum
        lt = size(sysCCGS2split{jsys}{isys}); %length of cellpair vector
        pre = zeros(lt(1),1);
        post = zeros(lt(1),1);
        
        for imat = 1: lt(1)
            pre(imat) = sum(sysCCGS2split{jsys}{isys}(imat,31:50));
            post(imat) = sum(sysCCGS2split{jsys}{isys}(imat,52:71));
        end
        
        if (normal)
            S2sessbias{isys,jsys} = (post - pre)./ (post + pre);
        else
            S2sessbias{isys,jsys} = (post - pre);
        end %if
    end % for jsys
end % for isys



