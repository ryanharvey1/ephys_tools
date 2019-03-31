%% Radial Arm Maze 
% Group 1 Probe Trial 11/6/15, ALL ARM CHOICES INCLUDED
% ST01
% w=[4 7 2];
w=[4 2];  % Getting rid of 3 and 7
w=createvector(w);
alpha_rad=armtoradians(w);

% ST02
% y=[4 6 2 4 7];
y=[4 6 2 4];
y=createvector(y);
beta_rad=armtoradians(y);

% ST04
% g=[5 7 4 2 4];
g=[5 4 2 4];
g=createvector(g);
gamma_rad=armtoradians(g);

% ST03
d=[5 2 8]; 
d=createvector(d);
delta_rad=armtoradians(d);

% ST07
% e=[3 8 5 2 1 7];
e=[8 5 2 1];
e=createvector(e);
epsilon_rad=armtoradians(e);

% ST08
% z=[7 2 8 4 2 6];
z=[2 8 4 2 6];
z=createvector(z);
zeta_rad=armtoradians(z);

% Normalized around 0, All choices
% ST01,02,04:[4 2 4 6 2 4 5 4 2 4 6] Normalized below in vector
all=[8 6 8 2 6 8 1 8 6 8 2 5 2 8 8 5 2 1 2 8 4 2 6];
all=createvector(all);
average_rad=armtoradians(all);

figure(1)
subplot(2,6,1)
circ_plot(alpha_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST01')
subplot(2,6,7)
circ_plot(alpha_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST01')
subplot(2,6,2)
circ_plot(beta_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST02')
subplot(2,6,8)
circ_plot(beta_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST02')
subplot(2,6,3)
circ_plot(gamma_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST04')
subplot(2,6,9)
circ_plot(gamma_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST04')
subplot(2,6,4)
circ_plot(delta_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST03')
subplot(2,6,10)
circ_plot(delta_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST03')
subplot(2,6,5)
circ_plot(epsilon_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST07')
subplot(2,6,11)
circ_plot(epsilon_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST07')
subplot(2,6,6)
circ_plot(zeta_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST08')
subplot(2,6,12)
circ_plot(zeta_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST08')

figure (2)
subplot(1,2,1)
circ_plot(average_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('Average')
subplot(1,2,2)
circ_plot(average_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('Average Hist')

figure (3)
compass(average_rad);

figure (4)
h = polar(average_rad,'.');
set(h,'markersize',12)

fprintf('\t\t\t\t\tALPHA\tBETA\tGAMMA\tDELTA\tEPSILON\tZETA\tAVERAGE\n')

alpha_bar = circ_mean(alpha_rad);
beta_bar = circ_mean(beta_rad);
gamma_bar=circ_mean(gamma_rad);
delta_bar=circ_mean(delta_rad);
epsilon_bar=circ_mean(epsilon_rad);
zeta_bar=circ_mean(zeta_rad);
average_bar=circ_mean(average_rad);

fprintf('Mean resultant vector:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar beta_bar gamma_bar delta_bar epsilon_bar zeta_bar average_bar]))

alpha_hat = circ_median(alpha_rad);
beta_hat = circ_median(beta_rad);
gamma_hat=circ_median(gamma_rad);
delta_hat=circ_median(delta_rad);
epsilon_hat=circ_median(epsilon_rad);
zeta_hat=circ_median(zeta_rad);
average_hat=circ_median(average_rad);

fprintf('Median:\t\t\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_hat beta_hat gamma_hat delta_hat epsilon_hat zeta_hat average_hat]))

R_alpha = circ_r(alpha_rad);
R_beta = circ_r(beta_rad);
R_gamma=circ_r(gamma_rad);
R_delta=circ_r(delta_rad);
R_epsilon=circ_r(epsilon_rad);
R_zeta=circ_r(zeta_rad);
R_average=circ_r(average_rad);

fprintf('R Length:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[R_alpha R_beta R_gamma R_delta R_epsilon R_zeta R_average])

S_alpha = circ_var(alpha_rad);
S_beta = circ_var(beta_rad);
S_gamma=circ_var(gamma_rad);
S_delta=circ_var(delta_rad);
S_epsilon=circ_var(epsilon_rad);
S_zeta=circ_var(zeta_rad);
S_average=circ_var(average_rad);

fprintf('Variance:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[S_alpha S_beta S_gamma S_delta S_epsilon S_zeta S_average])

[s_alpha s0_alpha] = circ_std(alpha_rad);
[s_beta s0_beta] = circ_std(beta_rad);
[s_gamma s0_gamma] = circ_std(gamma_rad);
[s_delta s0_delta] = circ_std(delta_rad);
[s_epsilon s0_epsilon] = circ_std(epsilon_rad);
[s_zeta s0_zeta] = circ_std(zeta_rad);
[s_average s0_average] = circ_std(average_rad);

fprintf('Standard deviation:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[s_alpha s_beta s_gamma s_delta s_epsilon s_zeta s_average])
fprintf('Standard deviation 0:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[s0_alpha s0_beta s0_gamma s0_delta s0_epsilon s0_zeta s0_average])

b_alpha = circ_skewness(alpha_rad);
b_beta = circ_skewness(beta_rad);
b_gamma = circ_skewness(gamma_rad);
b_delta = circ_skewness(delta_rad);
b_epsilon = circ_skewness(epsilon_rad);
b_zeta = circ_skewness(zeta_rad);
b_average = circ_skewness(average_rad);

fprintf('Skewness:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[b_alpha b_beta b_gamma b_delta b_epsilon b_zeta b_average])

k_alpha = circ_kurtosis(alpha_rad);
k_beta = circ_kurtosis(beta_rad);
k_gamma = circ_kurtosis(gamma_rad);
k_delta = circ_kurtosis(delta_rad);
k_epsilon = circ_kurtosis(epsilon_rad);
k_zeta = circ_kurtosis(zeta_rad);
k_average = circ_kurtosis(average_rad);

fprintf('Kurtosis:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[k_alpha k_beta k_gamma k_delta k_epsilon k_zeta k_average])

fprintf('\n\n')

fprintf('Inferential Statistics\n\nTests for Uniformity\n')
fprintf('\t\t\t\t\tALPHA\tBETA\tGAMMA\tDELTA\tEPSILON\tZETA\tAVERAGE\n')

% Rayleigh test
p_alpha = circ_rtest(alpha_rad);
p_beta = circ_rtest(beta_rad);
p_gamma = circ_rtest(gamma_rad);
p_delta = circ_rtest(delta_rad);
p_epsilon = circ_rtest(epsilon_rad);
p_zeta = circ_rtest(zeta_rad);
p_average = circ_rtest(average_rad);

fprintf('Rayleigh Test:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_epsilon p_zeta p_average])
             
% Omnibus test
o_alpha = circ_otest(alpha_rad);
o_beta = circ_otest(beta_rad);
o_gamma = circ_otest(gamma_rad);
o_delta = circ_otest(delta_rad);
o_epsilon = circ_otest(epsilon_rad);
o_zeta = circ_otest(zeta_rad);
o_average = circ_otest(average_rad);
fprintf('Omnibus Test:\t\t\t\t %.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[o_alpha o_beta o_gamma o_delta o_epsilon o_zeta o_average])

% Rao's spacing test
p_alpha = circ_raotest(alpha_rad);
p_beta = circ_raotest(beta_rad);
p_gamma = circ_raotest(gamma_rad);
p_delta = circ_raotest(delta_rad);
p_epsilon = circ_raotest(epsilon_rad);
p_zeta = circ_raotest(zeta_rad);
p_average = circ_raotest(average_rad);
fprintf('Rao Spacing Test:\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_epsilon p_zeta p_average])

% V test
p_alpha = circ_vtest(alpha_rad,circ_ang2rad(0));
p_beta = circ_vtest(beta_rad,circ_ang2rad(0));
p_gamma = circ_vtest(gamma_rad,circ_ang2rad(0));
p_delta = circ_vtest(delta_rad,circ_ang2rad(0));
p_epsilon = circ_vtest(epsilon_rad,circ_ang2rad(0));
p_zeta = circ_vtest(zeta_rad,circ_ang2rad(0));
p_average = circ_vtest(average_rad,circ_ang2rad(0));
fprintf('V Test (r = 0):\t\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_epsilon p_zeta p_average])

fprintf('\nTests concerning Mean and Median angle\n')
fprintf('\nDo not trust Cofidence means***********\n') % Had to create custom code for this
fprintf('\nThey are probably wrong\n')
% 95 percent confidence intervals for mean direction
w = ones(size(alpha_rad));
r = circ_r(alpha_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_alpha = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(beta_rad));
r = circ_r(beta_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_beta = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(gamma_rad));
r = circ_r(gamma_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_gamma = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(delta_rad));
r = circ_r(delta_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_delta = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(epsilon_rad));
r = circ_r(epsilon_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_epsilon = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(zeta_rad));
r = circ_r(zeta_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_zeta = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(average_rad));
r = circ_r(average_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_average = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

%t_alpha = circ_confmean(alpha_rad,0.05);
%t_beta = circ_confmean(beta_rad,0.05);
%t_gamma = circ_confmean(gamma_rad,0.05);
%t_delta = circ_confmean(delta_rad,0.05);
%t_epsilon = circ_confmean(epsilon_rad,0.05);
%t_zeta = circ_confmean(zeta_rad,0.05);
%t_average = circ_confmean(average_rad,0.05);

fprintf('Mean, up 95 perc. CI:\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar+t_alpha beta_bar+t_beta gamma_bar+t_gamma delta_bar+t_delta epsilon_bar+t_epsilon zeta_bar+t_zeta average_bar+t_average]))
fprintf('Mean, low 95 perc. CI:\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([2*pi+alpha_bar-t_alpha 2*pi+beta_bar-t_beta 2*pi+gamma_bar-t_gamma 2*pi+delta_bar-t_delta 2*pi+epsilon_bar-t_epsilon 2*pi+zeta_bar-t_zeta 2*pi+average_bar-t_average])) 

%% Radial Arm Maze 
% Group 1 Probe Trial ********JUST FIRST ARM CHOICE*******
%ST01
w=[4];
w(w==1) = 0;
w(w==2) = 45;
w(w==3) = 90;
w(w==4) = 135;
w(w==5) = 180;
w(w==6) = 225;
w(w==7) = 270;
w(w==8) = 315;
% ST02
y=[4];
y(y==1) = 0;
y(y==2) = 45;
y(y==3) = 90;
y(y==4) = 135;
y(y==5) = 180;
y(y==6) = 225;
y(y==7) = 270;
y(y==8) = 315;
% ST04
g=[5];
g(g==1) = 0;
g(g==2) = 45;
g(g==3) = 90;
g(g==4) = 135;
g(g==5) = 180;
g(g==6) = 225;
g(g==7) = 270;
g(g==8) = 315;
% ST03
d=[5]; 
d(d==1) = 0;
d(d==2) = 45;
d(d==3) = 90;
d(d==4) = 135;
d(d==5) = 180;
d(d==6) = 225;
d(d==7) = 270;
d(d==8) = 315;
% ST07
e=[3];
e(e==1) = 0;
e(e==2) = 45;
e(e==3) = 90;
e(e==4) = 135;
e(e==5) = 180;
e(e==6) = 225;
e(e==7) = 270;
e(e==8) = 315;
% ST08
z=[7];
z(z==1) = 0;
z(z==2) = 45;
z(z==3) = 90;
z(z==4) = 135;
z(z==5) = 180;
z(z==6) = 225;
z(z==7) = 270;
z(z==8) = 315;

% Normalized around 0, All choices
% ST01,02,04: Normalized below in vector
average=[8 8 1 5 3 7];
average(average==1) = 0;
average(average==2) = 45;
average(average==3) = 90;
average(average==4) = 135;
average(average==5) = 180;
average(average==6) = 225;
average(average==7) = 270;
average(average==8) = 315;

average_deg=average';
average_rad = circ_ang2rad(average_deg); 

alpha_deg=w';
alpha_rad = circ_ang2rad(alpha_deg); 

beta_deg =y';
beta_rad = circ_ang2rad(beta_deg); 

gamma_deg =g';
gamma_rad = circ_ang2rad(gamma_deg); 

delta_deg =d';
delta_rad = circ_ang2rad(delta_deg); 

epsilon_deg =e';
epsilon_rad = circ_ang2rad(epsilon_deg);

zeta_deg =z';
zeta_rad = circ_ang2rad(zeta_deg);


figure(3)
subplot(2,6,1)
circ_plot(alpha_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST01')
subplot(2,6,7)
circ_plot(alpha_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST01')
subplot(2,6,2)
circ_plot(beta_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST02')
subplot(2,6,8)
circ_plot(beta_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST02')
subplot(2,6,3)
circ_plot(gamma_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST04')
subplot(2,6,9)
circ_plot(gamma_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST04')
subplot(2,6,4)
circ_plot(delta_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST03')
subplot(2,6,10)
circ_plot(delta_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST03')
subplot(2,6,5)
circ_plot(epsilon_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST07')
subplot(2,6,11)
circ_plot(epsilon_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST07')
subplot(2,6,6)
circ_plot(zeta_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST08')
subplot(2,6,12)
circ_plot(zeta_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST08')

figure (4)
subplot(1,2,1)
circ_plot(average_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('Average')
subplot(1,2,2)
circ_plot(average_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('Average Hist')

fprintf('\t\t\t\t\tALPHA\tBETA\tGAMMA\tDELTA\tEPSILON\tZETA\tAVERAGE\n')

alpha_bar = circ_mean(alpha_rad);
beta_bar = circ_mean(beta_rad);
gamma_bar=circ_mean(gamma_rad);
delta_bar=circ_mean(delta_rad);
epsilon_bar=circ_mean(epsilon_rad);
zeta_bar=circ_mean(zeta_rad);
average_bar=circ_mean(average_rad);

fprintf('Mean resultant vector:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar beta_bar gamma_bar delta_bar epsilon_bar zeta_bar average_bar]))

alpha_hat = circ_median(alpha_rad);
beta_hat = circ_median(beta_rad);
gamma_hat=circ_median(gamma_rad);
delta_hat=circ_median(delta_rad);
epsilon_hat=circ_median(epsilon_rad);
zeta_hat=circ_median(zeta_rad);
average_hat=circ_median(average_rad);

fprintf('Median:\t\t\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_hat beta_hat gamma_hat delta_hat epsilon_hat zeta_hat average_hat]))

R_alpha = circ_r(alpha_rad);
R_beta = circ_r(beta_rad);
R_gamma=circ_r(gamma_rad);
R_delta=circ_r(delta_rad);
R_epsilon=circ_r(epsilon_rad);
R_zeta=circ_r(zeta_rad);
R_average=circ_r(average_rad);

fprintf('R Length:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[R_alpha R_beta R_gamma R_delta R_epsilon R_zeta R_average])

S_alpha = circ_var(alpha_rad);
S_beta = circ_var(beta_rad);
S_gamma=circ_var(gamma_rad);
S_delta=circ_var(delta_rad);
S_epsilon=circ_var(epsilon_rad);
S_zeta=circ_var(zeta_rad);
S_average=circ_var(average_rad);

fprintf('Variance:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[S_alpha S_beta S_gamma S_delta S_epsilon S_zeta S_average])

[s_alpha s0_alpha] = circ_std(alpha_rad);
[s_beta s0_beta] = circ_std(beta_rad);
[s_gamma s0_gamma] = circ_std(gamma_rad);
[s_delta s0_delta] = circ_std(delta_rad);
[s_epsilon s0_epsilon] = circ_std(epsilon_rad);
[s_zeta s0_zeta] = circ_std(zeta_rad);
[s_average s0_average] = circ_std(average_rad);

fprintf('Standard deviation:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[s_alpha s_beta s_gamma s_delta s_epsilon s_zeta s_average])
fprintf('Standard deviation 0:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[s0_alpha s0_beta s0_gamma s0_delta s0_epsilon s0_zeta s0_average])

b_alpha = circ_skewness(alpha_rad);
b_beta = circ_skewness(beta_rad);
b_gamma = circ_skewness(gamma_rad);
b_delta = circ_skewness(delta_rad);
b_epsilon = circ_skewness(epsilon_rad);
b_zeta = circ_skewness(zeta_rad);
b_average = circ_skewness(average_rad);

fprintf('Skewness:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[b_alpha b_beta b_gamma b_delta b_epsilon b_zeta b_average])

k_alpha = circ_kurtosis(alpha_rad);
k_beta = circ_kurtosis(beta_rad);
k_gamma = circ_kurtosis(gamma_rad);
k_delta = circ_kurtosis(delta_rad);
k_epsilon = circ_kurtosis(epsilon_rad);
k_zeta = circ_kurtosis(zeta_rad);
k_average = circ_kurtosis(average_rad);

fprintf('Kurtosis:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[k_alpha k_beta k_gamma k_delta k_epsilon k_zeta k_average])

fprintf('\n\n')

fprintf('Inferential Statistics\n\nTests for Uniformity\n')
fprintf('\t\t\t\t\tALPHA\tBETA\tGAMMA\tDELTA\tEPSILON\tZETA\tAVERAGE\n')

% Rayleigh test
p_alpha = circ_rtest(alpha_rad);
p_beta = circ_rtest(beta_rad);
p_gamma = circ_rtest(gamma_rad);
p_delta = circ_rtest(delta_rad);
p_epsilon = circ_rtest(epsilon_rad);
p_zeta = circ_rtest(zeta_rad);
p_average = circ_rtest(average_rad);

fprintf('Rayleigh Test:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_epsilon p_zeta p_average])
             
% Omnibus test
o_alpha = circ_otest(alpha_rad);
o_beta = circ_otest(beta_rad);
o_gamma = circ_otest(gamma_rad);
o_delta = circ_otest(delta_rad);
o_epsilon = circ_otest(epsilon_rad);
o_zeta = circ_otest(zeta_rad);
o_average = circ_otest(average_rad);
fprintf('Omnibus Test:\t\t\t\t %.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[o_alpha o_beta o_gamma o_delta o_epsilon o_zeta o_average])

% Rao's spacing test
p_alpha = circ_raotest(alpha_rad);
p_beta = circ_raotest(beta_rad);
p_gamma = circ_raotest(gamma_rad);
p_delta = circ_raotest(delta_rad);
p_epsilon = circ_raotest(epsilon_rad);
p_zeta = circ_raotest(zeta_rad);
p_average = circ_raotest(average_rad);
fprintf('Rao Spacing Test:\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_epsilon p_zeta p_average])

% V test
p_alpha = circ_vtest(alpha_rad,circ_ang2rad(0));
p_beta = circ_vtest(beta_rad,circ_ang2rad(0));
p_gamma = circ_vtest(gamma_rad,circ_ang2rad(0));
p_delta = circ_vtest(delta_rad,circ_ang2rad(0));
p_epsilon = circ_vtest(epsilon_rad,circ_ang2rad(0));
p_zeta = circ_vtest(zeta_rad,circ_ang2rad(0));
p_average = circ_vtest(average_rad,circ_ang2rad(0));
fprintf('V Test (r = 0):\t\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_epsilon p_zeta p_average])

fprintf('\nTests concerning Mean and Median angle\n')
fprintf('\nDo not trust Cofidence means***********\n') % Had to create custom code for this
fprintf('\nThey are probably wrong\n')
% 95 percent confidence intervals for mean direction
w = ones(size(alpha_rad));
r = circ_r(alpha_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_alpha = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(beta_rad));
r = circ_r(beta_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_beta = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(gamma_rad));
r = circ_r(gamma_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_gamma = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(delta_rad));
r = circ_r(delta_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_delta = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(epsilon_rad));
r = circ_r(epsilon_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_epsilon = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(zeta_rad));
r = circ_r(zeta_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_zeta = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(average_rad));
r = circ_r(average_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_average = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

%t_alpha = circ_confmean(alpha_rad,0.05);
%t_beta = circ_confmean(beta_rad,0.05);
%t_gamma = circ_confmean(gamma_rad,0.05);
%t_delta = circ_confmean(delta_rad,0.05);
%t_epsilon = circ_confmean(epsilon_rad,0.05);
%t_zeta = circ_confmean(zeta_rad,0.05);
%t_average = circ_confmean(average_rad,0.05);

fprintf('Mean, up 95 perc. CI:\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar+t_alpha beta_bar+t_beta gamma_bar+t_gamma delta_bar+t_delta epsilon_bar+t_epsilon zeta_bar+t_zeta average_bar+t_average]))
fprintf('Mean, low 95 perc. CI:\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([2*pi+alpha_bar-t_alpha 2*pi+beta_bar-t_beta 2*pi+gamma_bar-t_gamma 2*pi+delta_bar-t_delta 2*pi+epsilon_bar-t_epsilon 2*pi+zeta_bar-t_zeta 2*pi+average_bar-t_average])) 

%% Radial Arm Maze 
% Group 1 Probe Trial ********JUST FIRST TWO ARM CHOICES*******
%ST01
w=[4 7];
w(w==1) = 0;
w(w==2) = 45;
w(w==3) = 90;
w(w==4) = 135;
w(w==5) = 180;
w(w==6) = 225;
w(w==7) = 270;
w(w==8) = 315;
% ST02
y=[4 6];
y(y==1) = 0;
y(y==2) = 45;
y(y==3) = 90;
y(y==4) = 135;
y(y==5) = 180;
y(y==6) = 225;
y(y==7) = 270;
y(y==8) = 315;
% ST04
g=[5 7];
g(g==1) = 0;
g(g==2) = 45;
g(g==3) = 90;
g(g==4) = 135;
g(g==5) = 180;
g(g==6) = 225;
g(g==7) = 270;
g(g==8) = 315;
% ST03
d=[5 2]; 
d(d==1) = 0;
d(d==2) = 45;
d(d==3) = 90;
d(d==4) = 135;
d(d==5) = 180;
d(d==6) = 225;
d(d==7) = 270;
d(d==8) = 315;
% ST07
e=[3 8];
e(e==1) = 0;
e(e==2) = 45;
e(e==3) = 90;
e(e==4) = 135;
e(e==5) = 180;
e(e==6) = 225;
e(e==7) = 270;
e(e==8) = 315;
% ST08
z=[7 2];
z(z==1) = 0;
z(z==2) = 45;
z(z==3) = 90;
z(z==4) = 135;
z(z==5) = 180;
z(z==6) = 225;
z(z==7) = 270;
z(z==8) = 315;

% Normalized around 0, All choices
% ST01,02,04: Normalized below in vector
average=[8 3 8 2 1 3 5 2 3 8 7 2];
average(average==1) = 0;
average(average==2) = 45;
average(average==3) = 90;
average(average==4) = 135;
average(average==5) = 180;
average(average==6) = 225;
average(average==7) = 270;
average(average==8) = 315;

average_deg=average';
average_rad = circ_ang2rad(average_deg); 

alpha_deg=w';
alpha_rad = circ_ang2rad(alpha_deg); 

beta_deg =y';
beta_rad = circ_ang2rad(beta_deg); 

gamma_deg =g';
gamma_rad = circ_ang2rad(gamma_deg); 

delta_deg =d';
delta_rad = circ_ang2rad(delta_deg); 

epsilon_deg =e';
epsilon_rad = circ_ang2rad(epsilon_deg);

zeta_deg =z';
zeta_rad = circ_ang2rad(zeta_deg);


figure(5)
subplot(2,6,1)
circ_plot(alpha_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST01')
subplot(2,6,7)
circ_plot(alpha_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST01')
subplot(2,6,2)
circ_plot(beta_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST02')
subplot(2,6,8)
circ_plot(beta_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST02')
subplot(2,6,3)
circ_plot(gamma_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST04')
subplot(2,6,9)
circ_plot(gamma_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST04')
subplot(2,6,4)
circ_plot(delta_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST03')
subplot(2,6,10)
circ_plot(delta_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST03')
subplot(2,6,5)
circ_plot(epsilon_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST07')
subplot(2,6,11)
circ_plot(epsilon_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST07')
subplot(2,6,6)
circ_plot(zeta_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('ST08')
subplot(2,6,12)
circ_plot(zeta_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('ST08')

figure (6)
subplot(1,2,1)
circ_plot(average_rad,'pretty','bo',true,'linewidth',2,'color','r'),
title('Average')
subplot(1,2,2)
circ_plot(average_rad,'hist',[],20,true,true,'linewidth',2,'color','r')
title('Average Hist')

fprintf('\t\t\t\t\tALPHA\tBETA\tGAMMA\tDELTA\tEPSILON\tZETA\tAVERAGE\n')

alpha_bar = circ_mean(alpha_rad);
beta_bar = circ_mean(beta_rad);
gamma_bar=circ_mean(gamma_rad);
delta_bar=circ_mean(delta_rad);
epsilon_bar=circ_mean(epsilon_rad);
zeta_bar=circ_mean(zeta_rad);
average_bar=circ_mean(average_rad);

fprintf('Mean resultant vector:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar beta_bar gamma_bar delta_bar epsilon_bar zeta_bar average_bar]))

alpha_hat = circ_median(alpha_rad);
beta_hat = circ_median(beta_rad);
gamma_hat=circ_median(gamma_rad);
delta_hat=circ_median(delta_rad);
epsilon_hat=circ_median(epsilon_rad);
zeta_hat=circ_median(zeta_rad);
average_hat=circ_median(average_rad);

fprintf('Median:\t\t\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_hat beta_hat gamma_hat delta_hat epsilon_hat zeta_hat average_hat]))

R_alpha = circ_r(alpha_rad);
R_beta = circ_r(beta_rad);
R_gamma=circ_r(gamma_rad);
R_delta=circ_r(delta_rad);
R_epsilon=circ_r(epsilon_rad);
R_zeta=circ_r(zeta_rad);
R_average=circ_r(average_rad);

fprintf('R Length:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[R_alpha R_beta R_gamma R_delta R_epsilon R_zeta R_average])

S_alpha = circ_var(alpha_rad);
S_beta = circ_var(beta_rad);
S_gamma=circ_var(gamma_rad);
S_delta=circ_var(delta_rad);
S_epsilon=circ_var(epsilon_rad);
S_zeta=circ_var(zeta_rad);
S_average=circ_var(average_rad);

fprintf('Variance:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[S_alpha S_beta S_gamma S_delta S_epsilon S_zeta S_average])

[s_alpha s0_alpha] = circ_std(alpha_rad);
[s_beta s0_beta] = circ_std(beta_rad);
[s_gamma s0_gamma] = circ_std(gamma_rad);
[s_delta s0_delta] = circ_std(delta_rad);
[s_epsilon s0_epsilon] = circ_std(epsilon_rad);
[s_zeta s0_zeta] = circ_std(zeta_rad);
[s_average s0_average] = circ_std(average_rad);

fprintf('Standard deviation:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[s_alpha s_beta s_gamma s_delta s_epsilon s_zeta s_average])
fprintf('Standard deviation 0:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[s0_alpha s0_beta s0_gamma s0_delta s0_epsilon s0_zeta s0_average])

b_alpha = circ_skewness(alpha_rad);
b_beta = circ_skewness(beta_rad);
b_gamma = circ_skewness(gamma_rad);
b_delta = circ_skewness(delta_rad);
b_epsilon = circ_skewness(epsilon_rad);
b_zeta = circ_skewness(zeta_rad);
b_average = circ_skewness(average_rad);

fprintf('Skewness:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[b_alpha b_beta b_gamma b_delta b_epsilon b_zeta b_average])

k_alpha = circ_kurtosis(alpha_rad);
k_beta = circ_kurtosis(beta_rad);
k_gamma = circ_kurtosis(gamma_rad);
k_delta = circ_kurtosis(delta_rad);
k_epsilon = circ_kurtosis(epsilon_rad);
k_zeta = circ_kurtosis(zeta_rad);
k_average = circ_kurtosis(average_rad);

fprintf('Kurtosis:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[k_alpha k_beta k_gamma k_delta k_epsilon k_zeta k_average])

fprintf('\n\n')

fprintf('Inferential Statistics\n\nTests for Uniformity\n')
fprintf('\t\t\t\t\tALPHA\tBETA\tGAMMA\tDELTA\tEPSILON\tZETA\tAVERAGE\n')

% Rayleigh test
p_alpha = circ_rtest(alpha_rad);
p_beta = circ_rtest(beta_rad);
p_gamma = circ_rtest(gamma_rad);
p_delta = circ_rtest(delta_rad);
p_epsilon = circ_rtest(epsilon_rad);
p_zeta = circ_rtest(zeta_rad);
p_average = circ_rtest(average_rad);

fprintf('Rayleigh Test:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_epsilon p_zeta p_average])
             
% Omnibus test
o_alpha = circ_otest(alpha_rad);
o_beta = circ_otest(beta_rad);
o_gamma = circ_otest(gamma_rad);
o_delta = circ_otest(delta_rad);
o_epsilon = circ_otest(epsilon_rad);
o_zeta = circ_otest(zeta_rad);
o_average = circ_otest(average_rad);
fprintf('Omnibus Test:\t\t\t\t %.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[o_alpha o_beta o_gamma o_delta o_epsilon o_zeta o_average])

% Rao's spacing test
p_alpha = circ_raotest(alpha_rad);
p_beta = circ_raotest(beta_rad);
p_gamma = circ_raotest(gamma_rad);
p_delta = circ_raotest(delta_rad);
p_epsilon = circ_raotest(epsilon_rad);
p_zeta = circ_raotest(zeta_rad);
p_average = circ_raotest(average_rad);
fprintf('Rao Spacing Test:\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_epsilon p_zeta p_average])

% V test
p_alpha = circ_vtest(alpha_rad,circ_ang2rad(0));
p_beta = circ_vtest(beta_rad,circ_ang2rad(0));
p_gamma = circ_vtest(gamma_rad,circ_ang2rad(0));
p_delta = circ_vtest(delta_rad,circ_ang2rad(0));
p_epsilon = circ_vtest(epsilon_rad,circ_ang2rad(0));
p_zeta = circ_vtest(zeta_rad,circ_ang2rad(0));
p_average = circ_vtest(average_rad,circ_ang2rad(0));
fprintf('V Test (r = 0):\t\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_epsilon p_zeta p_average])

fprintf('\nTests concerning Mean and Median angle\n')
fprintf('\nDo not trust Cofidence means***********\n') % Had to create custom code for this
fprintf('\nThey are probably wrong\n')
% 95 percent confidence intervals for mean direction
w = ones(size(alpha_rad));
r = circ_r(alpha_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_alpha = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(beta_rad));
r = circ_r(beta_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_beta = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(gamma_rad));
r = circ_r(gamma_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_gamma = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(delta_rad));
r = circ_r(delta_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_delta = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(epsilon_rad));
r = circ_r(epsilon_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_epsilon = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(zeta_rad));
r = circ_r(zeta_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_zeta = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

w = ones(size(average_rad));
r = circ_r(average_rad,w,0,1);
xi = 0.05;
c2 = chi2inv((1-xi),1);
n = sum(w,1);
R = n.*r;
t = zeros(size(r));
i = 1:numel(r);
sqrt(c2/2/n(i));
t_average = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));

%t_alpha = circ_confmean(alpha_rad,0.05);
%t_beta = circ_confmean(beta_rad,0.05);
%t_gamma = circ_confmean(gamma_rad,0.05);
%t_delta = circ_confmean(delta_rad,0.05);
%t_epsilon = circ_confmean(epsilon_rad,0.05);
%t_zeta = circ_confmean(zeta_rad,0.05);
%t_average = circ_confmean(average_rad,0.05);

fprintf('Mean, up 95 perc. CI:\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar+t_alpha beta_bar+t_beta gamma_bar+t_gamma delta_bar+t_delta epsilon_bar+t_epsilon zeta_bar+t_zeta average_bar+t_average]))
fprintf('Mean, low 95 perc. CI:\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([2*pi+alpha_bar-t_alpha 2*pi+beta_bar-t_beta 2*pi+gamma_bar-t_gamma 2*pi+delta_bar-t_delta 2*pi+epsilon_bar-t_epsilon 2*pi+zeta_bar-t_zeta 2*pi+average_bar-t_average])) 
