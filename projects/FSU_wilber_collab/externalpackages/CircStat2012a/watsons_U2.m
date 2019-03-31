function U2=watsons_U2(A1,A2)
%
% U2=watsons_U2(A1,A2)
%
% Computes Watson's U2 statistic for nonparametric 2-sample testing of
% circular data, accommodating ties. The code is derived from eq. 27.17
% Zar (1999) and its performance was verified using the numerical examples
% 27.10 and 27.11 from that reference.
%
% Inputs:
% A1, A2:   vectors containing angles (in degrees or radians, the unit does
%           not matter)
% 
% Outputs:
% U2:       Watson's U2 statistic
%
% Significance tables for U2 have been published, e.g. in Zar (1999).
% Alternatively, an ad hoc permutation test can be used.
%
% References:
% Zar JH. Biostatistical Analysis. 4th ed., 1999.
% Chapter 27: Circular Distributions: Hypothesis Testing.
% Upper Saddle River, NJ: Prentice Hall.
%
% pierre.megevand@gmail.com
%
% See also WATSONS_U2_PERM_TEST

[a1,t1,m1,n1,m1_n1]=uniques_ties_cumuls_relfreqs(A1);
[a2,t2,m2,n2,m2_n2]=uniques_ties_cumuls_relfreqs(A2);

n=n1+n2;

k=numel(unique([A1;A2]));
table=[unique([A1;A2]) zeros(k,3)];
for i=1:size(table,1)
    [tf1,loc1]=ismember(table(i,1),a1);
    [tf2,loc2]=ismember(table(i,1),a2);
    if tf1
        table(i,2)=table(i,2)+m1_n1(loc1);
        table(i,4)=table(i,4)+t1(loc1);
    else
        if i>1
            table(i,2)=table(i-1,2);
        end
    end
    if tf2
        table(i,3)=table(i,3)+m2_n2(loc2);
        table(i,4)=table(i,4)+t2(loc2);
    else
        if i>1
            table(i,3)=table(i-1,3);
        end
    end
end

d=table(:,2)-table(:,3);
t=table(:,4);
td=sum(t.*d);
td2=sum(t.*d.^2);

U2=((n1.*n2)./(n.^2)).*(td2-(td.^2./n));

    function [a,t,m,n,m_n]=uniques_ties_cumuls_relfreqs(A)
        a=unique(A);
        t=zeros(size(a));
        m=zeros(size(a));
        for j=1:length(a)
            t(j)=numel(A(A==a(j)));
            m(j)=sum(t(1:j));
        end
        n=m(end);
        m_n=m./n;
    end

end

