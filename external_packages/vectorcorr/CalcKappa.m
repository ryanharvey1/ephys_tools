function  Result = CalcKappa(IntVal, flag)

%CalcKappa.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

%Look up Kappa values from R-bar values if flag = 0
%Look up R-bar values from Kappa values if flag = 1

% Variables input:
%   IntVal: value to be interpolated
%   flag:   0=interpolate Kappa from R-bar   1=interpolate R-bar
% Output variable:
%   Result: interpolated value

%TableK contains R-bar in col.1, KappaHat in col.2
%  (from Mardia, 1972, p. 298; and Fisher, 1993, p.224)

TableK = [ 0       0; .01 .02;    .02 .04001; .03 .06003; .04 .08006;
          .05 .10013; .06 .12022; .07 .14034; .08 .16051; .09 .18073;
          .10 .20101; .11 .22134; .12 .24175; .13 .26223; .14 .28279; 
          .15 .30344; .16 .32419; .17 .34503; .18 .36599; .19 .38707;
          .20 .40828; .21 .42962; .22 .45110; .23 .47273; .24 .49453;
          .25 .51649; .26 .53863; .27 .56097; .28 .58350; .29 .60625;
          .30 .62922; .31 .65242; .32 .67587; .33 .69958; .34 .72356;
          .35 .74783; .36 .77241; .37 .79730; .38 .82253; .39 .84812;
          .40 .87408; .41 .90043; .42 .92720; .43 .95440; .44 .98207;
          .45 1.0102; .46 1.0389; .47 1.0681; .48 1.0979; .49 1.1283;
          .50 1.1593; .51 1.1911; .52 1.2235; .53 1.2567; .54 1.2907;
          .55 1.3257; .56 1.3616; .57 1.3984; .58 1.4364; .59 1.4754;
          .60 1.5157; .61 1.5574; .62 1.6004; .63 1.6451; .64 1.6913;
          .65 1.7395; .66 1.7895; .67 1.8418; .68 1.8964; .69 1.9536;
          .70 2.0136; .71 2.0769; .72 2.1436; .73 2.2142; .74 2.2893; 
          .75 2.3693; .76 2.4549; .77 2.5469; .78 2.6461; .79 2.7538;
          .80 2.8713; .81 3.0002; .82 3.1426; .83 3.3011; .84 3.4790;
          .85 3.6804; .86 3.9107; .87 4.1770; .88 4.4888; .89 4.8587;
          .90 5.3047; .91 5.8522; .92 6.5394; .93 7.4257; .94 8.6104;
          .95 10.272; .96 12.766; .97 16.927; .98 25.252; .99 50.242;
          .995 100;   .999  500.; 1.00  5000];

if flag == 0
    Result = interp1(TableK(:,1), TableK(:,2), IntVal, 'cubic');
else
    Result = interp1(TableK(:,2), TableK(:,1), IntVal, 'cubic');
end

      