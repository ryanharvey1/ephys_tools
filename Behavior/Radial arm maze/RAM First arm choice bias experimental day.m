% Radial Arm Maze First arm choice bias experimental day
% Group 4/5 Muscimol test 

% Control

% Convert arm choices into degrees 
% Trial 1
w=[5 1 1 5 2 5 4 4 4 4 2 7];
w(w==1) = 0;
w(w==2) = 45;
w(w==3) = 90;
w(w==4) = 135;
w(w==5) = 180;
w(w==6) = 225;
w(w==7) = 270;
w(w==8) = 315;
% Trial 2
y=[1 1 1 1 8 1 1 8 1 8 2 1];
y(y==1) = 0;
y(y==2) = 45;
y(y==3) = 90;
y(y==4) = 135;
y(y==5) = 180;
y(y==6) = 225;
y(y==7) = 270;
y(y==8) = 315;
% Trial 3
g=[3 3 4 3 3 3 2 3 4 3 3 5];
g(g==1) = 0;
g(g==2) = 45;
g(g==3) = 90;
g(g==4) = 135;
g(g==5) = 180;
g(g==6) = 225;
g(g==7) = 270;
g(g==8) = 315;
% Trial 4
d=[7 7 4 7 7 1 7 3 7 6 7 7]; 
d(d==1) = 0;
d(d==2) = 45;
d(d==3) = 90;
d(d==4) = 135;
d(d==5) = 180;
d(d==6) = 225;
d(d==7) = 270;
d(d==8) = 315;

% total arm choices 
% averagematrix=[];
average=[5 1 1 5 2 5 4 4 4 4 2 7 1	1	1	1	8	1	1	8	1	8	2	1 3	3	4	3	3	3	2	3	4	3	3	5 7	7	4	7	7	1	7	3	7	6	7	7];
average(average==1) = 0;
average(average==2) = 45;
average(average==3) = 90;
average(average==4) = 135;
average(average==5) = 180;
average(average==6) = 225;
average(average==7) = 270;
average(average==8) = 315;

% Now convert degrees into radians 
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

% Plot the circular directional heading 
figure(1)
subplot(2,4,1)
circ_plot(alpha_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Trial 1')
subplot(2,4,5)
circ_plot(alpha_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Trial 1')
subplot(2,4,2)
circ_plot(beta_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Trial 2')
subplot(2,4,6)
circ_plot(beta_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Trial 2')
subplot(2,4,3)
circ_plot(gamma_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Trial 3')
subplot(2,4,7)
circ_plot(gamma_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Trial 3')
subplot(2,4,4)
circ_plot(delta_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Trial 4')
subplot(2,4,8)
circ_plot(delta_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Trial 4')

figure (2);
subplot(1,2,1)
circ_plot(average_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Average')
subplot(1,2,2)
circ_plot(average_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Average Hist')

fprintf('\t\t\t\t\tALPHA\tBETA\tGAMMA\tDELTA\tAVERAGE\n')

alpha_bar = circ_mean(alpha_rad);
beta_bar = circ_mean(beta_rad);
gamma_bar=circ_mean(gamma_rad);
delta_bar=circ_mean(delta_rad);
average_bar=circ_mean(average_rad);

fprintf('Mean resultant vector:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar beta_bar gamma_bar delta_bar average_bar]))

alpha_hat = circ_median(alpha_rad);
beta_hat = circ_median(beta_rad);
gamma_hat=circ_median(gamma_rad);
delta_hat=circ_median(delta_rad);
average_hat=circ_median(average_rad);

fprintf('Median:\t\t\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_hat beta_hat gamma_hat delta_hat average_hat]))

R_alpha = circ_r(alpha_rad);
R_beta = circ_r(beta_rad);
R_gamma=circ_r(gamma_rad);
R_delta=circ_r(delta_rad);
R_average=circ_r(average_rad);

fprintf('R Length:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[R_alpha R_beta R_gamma R_delta R_average])

S_alpha = circ_var(alpha_rad);
S_beta = circ_var(beta_rad);
S_gamma=circ_var(gamma_rad);
S_delta=circ_var(delta_rad);
S_average=circ_var(average_rad);

fprintf('Variance:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[S_alpha S_beta S_gamma S_delta S_average])

[s_alpha s0_alpha] = circ_std(alpha_rad);
[s_beta s0_beta] = circ_std(beta_rad);
[s_gamma s0_gamma] = circ_std(gamma_rad);
[s_delta s0_delta] = circ_std(delta_rad);
[s_average s0_average] = circ_std(average_rad);

fprintf('Standard deviation:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[s_alpha s_beta s_gamma s_delta s_average])
fprintf('Standard deviation 0:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[s0_alpha s0_beta s0_gamma s0_delta s0_average])

b_alpha = circ_skewness(alpha_rad);
b_beta = circ_skewness(beta_rad);
b_gamma = circ_skewness(gamma_rad);
b_delta = circ_skewness(delta_rad);
b_average = circ_skewness(average_rad);

fprintf('Skewness:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[b_alpha b_beta b_gamma b_delta b_average])

k_alpha = circ_kurtosis(alpha_rad);
k_beta = circ_kurtosis(beta_rad);
k_gamma = circ_kurtosis(gamma_rad);
k_delta = circ_kurtosis(delta_rad);
k_average = circ_kurtosis(average_rad);

fprintf('Kurtosis:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[k_alpha k_beta k_gamma k_delta k_average])

fprintf('\n\n')

fprintf('Inferential Statistics\n\nTests for Uniformity\n')
fprintf('\t\t\t\t\tALPHA\tBETA\tGAMMA\tDELTA\tAVERAGE\n')

% Rayleigh test
p_alpha = circ_rtest(alpha_rad);
p_beta = circ_rtest(beta_rad);
p_gamma = circ_rtest(gamma_rad);
p_delta = circ_rtest(delta_rad);
p_average = circ_rtest(average_rad);

fprintf('Rayleigh Test:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_average])
             
% Omnibus test
o_alpha = circ_otest(alpha_rad);
o_beta = circ_otest(beta_rad);
o_gamma = circ_otest(gamma_rad);
o_delta = circ_otest(delta_rad);
o_average = circ_otest(average_rad);
fprintf('Omnibus Test:\t\t\t\t %.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[o_alpha o_beta o_gamma o_delta o_average])

% Rao's spacing test
p_alpha = circ_raotest(alpha_rad);
p_beta = circ_raotest(beta_rad);
p_gamma = circ_raotest(gamma_rad);
p_delta = circ_raotest(delta_rad);
p_average = circ_raotest(average_rad);
fprintf('Rao Spacing Test:\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_average])

% V test
p_alpha = circ_vtest(alpha_rad,circ_ang2rad(0));
p_beta = circ_vtest(beta_rad,circ_ang2rad(0));
p_gamma = circ_vtest(gamma_rad,circ_ang2rad(0));
p_delta = circ_vtest(delta_rad,circ_ang2rad(0));
p_average = circ_vtest(average_rad,circ_ang2rad(0));
fprintf('V Test (r = 0):\t\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_average])

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

fprintf('Mean, up 95 perc. CI:\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar+t_alpha beta_bar+t_beta gamma_bar+t_gamma delta_bar+t_delta average_bar+t_average]))
fprintf('Mean, low 95 perc. CI:\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([2*pi+alpha_bar-t_alpha 2*pi+beta_bar-t_beta 2*pi+gamma_bar-t_gamma 2*pi+delta_bar-t_delta 2*pi+average_bar-t_average])) 

%% ***************************************************************************************************
%*****************************************************************************************************
% Radial Arm Maze Firt arm choice bias

% Muscimol

% Convert arm choices into degrees 
% Trial 1
w=[1	5	6	6	5	5	4	5	8	4	6	4];
w(w==1) = 0;
w(w==2) = 45;
w(w==3) = 90;
w(w==4) = 135;
w(w==5) = 180;
w(w==6) = 225;
w(w==7) = 270;
w(w==8) = 315;
% Trial 2
y=[1	8	1	1	8	6	1	1	8	2		3];
y(y==1) = 0;
y(y==2) = 45;
y(y==3) = 90;
y(y==4) = 135;
y(y==5) = 180;
y(y==6) = 225;
y(y==7) = 270;
y(y==8) = 315;
% Trial 3
g=[4	3	7	4	2	2	2	3	2	4	4	3];
g(g==1) = 0;
g(g==2) = 45;
g(g==3) = 90;
g(g==4) = 135;
g(g==5) = 180;
g(g==6) = 225;
g(g==7) = 270;
g(g==8) = 315;
% Trial 4
d=[	1	1	8	7	6	6	7	8	6		5];
d(d==1) = 0;
d(d==2) = 45;
d(d==3) = 90;
d(d==4) = 135;
d(d==5) = 180;
d(d==6) = 225;
d(d==7) = 270;
d(d==8) = 315;

% total arm choices 
% averagematrix=(data([120:125],[1:4]));
average=[1	5	6	6	5	5	4	5	8	4	6	4 1	8	1	1	8	6	1	1	8	2		3 4	3	7	4	2	2	2	3	2	4	4	3 	1	1	8	7	6	6	7	8	6		5];
average(average==1) = 0;
average(average==2) = 45;
average(average==3) = 90;
average(average==4) = 135;
average(average==5) = 180;
average(average==6) = 225;
average(average==7) = 270;
average(average==8) = 315;

% Now convert degrees into radians 
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

% Plot the circular directional heading 
figure(3)
title('Muscimol')
subplot(2,4,1)
circ_plot(alpha_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Trial 1')
subplot(2,4,5)
circ_plot(alpha_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Trial 1')
subplot(2,4,2)
circ_plot(beta_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Trial 2')
subplot(2,4,6)
circ_plot(beta_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Trial 2')
subplot(2,4,3)
circ_plot(gamma_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Trial 3')
subplot(2,4,7)
circ_plot(gamma_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Trial 3')
subplot(2,4,4)
circ_plot(delta_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Trial 4')
subplot(2,4,8)
circ_plot(delta_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Trial 4')

figure (4);
subplot(1,2,1)
circ_plot(average_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Average')
subplot(1,2,2)
circ_plot(average_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Average Hist')

fprintf('\t\t\t\t\tALPHA\tBETA\tGAMMA\tDELTA\tAVERAGE\n')

alpha_bar = circ_mean(alpha_rad);
beta_bar = circ_mean(beta_rad);
gamma_bar=circ_mean(gamma_rad);
delta_bar=circ_mean(delta_rad);
average_bar=circ_mean(average_rad);

fprintf('Mean resultant vector:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar beta_bar gamma_bar delta_bar average_bar]))

alpha_hat = circ_median(alpha_rad);
beta_hat = circ_median(beta_rad);
gamma_hat=circ_median(gamma_rad);
delta_hat=circ_median(delta_rad);
average_hat=circ_median(average_rad);

fprintf('Median:\t\t\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_hat beta_hat gamma_hat delta_hat average_hat]))

R_alpha = circ_r(alpha_rad);
R_beta = circ_r(beta_rad);
R_gamma=circ_r(gamma_rad);
R_delta=circ_r(delta_rad);
R_average=circ_r(average_rad);

fprintf('R Length:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[R_alpha R_beta R_gamma R_delta R_average])

S_alpha = circ_var(alpha_rad);
S_beta = circ_var(beta_rad);
S_gamma=circ_var(gamma_rad);
S_delta=circ_var(delta_rad);
S_average=circ_var(average_rad);

fprintf('Variance:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[S_alpha S_beta S_gamma S_delta S_average])

[s_alpha s0_alpha] = circ_std(alpha_rad);
[s_beta s0_beta] = circ_std(beta_rad);
[s_gamma s0_gamma] = circ_std(gamma_rad);
[s_delta s0_delta] = circ_std(delta_rad);
[s_average s0_average] = circ_std(average_rad);

fprintf('Standard deviation:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[s_alpha s_beta s_gamma s_delta s_average])
fprintf('Standard deviation 0:\t\t\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[s0_alpha s0_beta s0_gamma s0_delta s0_average])

b_alpha = circ_skewness(alpha_rad);
b_beta = circ_skewness(beta_rad);
b_gamma = circ_skewness(gamma_rad);
b_delta = circ_skewness(delta_rad);
b_average = circ_skewness(average_rad);

fprintf('Skewness:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[b_alpha b_beta b_gamma b_delta b_average])

k_alpha = circ_kurtosis(alpha_rad);
k_beta = circ_kurtosis(beta_rad);
k_gamma = circ_kurtosis(gamma_rad);
k_delta = circ_kurtosis(delta_rad);
k_average = circ_kurtosis(average_rad);

fprintf('Kurtosis:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[k_alpha k_beta k_gamma k_delta k_average])

fprintf('\n\n')

fprintf('Inferential Statistics\n\nTests for Uniformity\n')
fprintf('\t\t\t\t\tALPHA\tBETA\tGAMMA\tDELTA\tAVERAGE\n')

% Rayleigh test
p_alpha = circ_rtest(alpha_rad);
p_beta = circ_rtest(beta_rad);
p_gamma = circ_rtest(gamma_rad);
p_delta = circ_rtest(delta_rad);
p_average = circ_rtest(average_rad);

fprintf('Rayleigh Test:\t\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_average])
             
% Omnibus test
o_alpha = circ_otest(alpha_rad);
o_beta = circ_otest(beta_rad);
o_gamma = circ_otest(gamma_rad);
o_delta = circ_otest(delta_rad);
o_average = circ_otest(average_rad);
fprintf('Omnibus Test:\t\t\t\t %.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[o_alpha o_beta o_gamma o_delta o_average])

% Rao's spacing test
p_alpha = circ_raotest(alpha_rad);
p_beta = circ_raotest(beta_rad);
p_gamma = circ_raotest(gamma_rad);
p_delta = circ_raotest(delta_rad);
p_average = circ_raotest(average_rad);
fprintf('Rao Spacing Test:\t\t\t%.2f \t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_average])

% V test
p_alpha = circ_vtest(alpha_rad,circ_ang2rad(0));
p_beta = circ_vtest(beta_rad,circ_ang2rad(0));
p_gamma = circ_vtest(gamma_rad,circ_ang2rad(0));
p_delta = circ_vtest(delta_rad,circ_ang2rad(0));
p_average = circ_vtest(average_rad,circ_ang2rad(0));
fprintf('V Test (r = 0):\t\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',[p_alpha p_beta p_gamma p_delta p_average])

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

fprintf('Mean, up 95 perc. CI:\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar+t_alpha beta_bar+t_beta gamma_bar+t_gamma delta_bar+t_delta average_bar+t_average]))
fprintf('Mean, low 95 perc. CI:\t\t\t %.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', circ_rad2ang([2*pi+alpha_bar-t_alpha 2*pi+beta_bar-t_beta 2*pi+gamma_bar-t_gamma 2*pi+delta_bar-t_delta 2*pi+average_bar-t_average])) 
