%% RAM TEST 4/5 all arm choices
% split up by experimental group

% import data
% data = xlsread('/Users/RyanHarvey/Dropbox/Lab/Projects/Anterior Thalamic inactivation/Data Sheets/First choice bias.xlsx');
% Convert arm choices into degrees 
% pooled and normalized around 1,4 
% Muscimol
w=[5 6	4	2	4	1	6	4	1	2	4	1	6	1	4	8	2	5	1	4	5	2	1	4	6	1	4	2	4	7	1	4	1	8	4	8	4	8	4	8	4	8	4	8	4	7	4	8	4	6	1	2	6	3	5	2	4	8	4	8	4	8	4	8	4	1	6	4	6	1	4	8	4	8	5	1	4	8	4	1	6	2	7	3	7	3	1	4	4	2	5	1	1	5	1	5	1	5	1	5	1	5	1	5	2	6	2	6	2	6	2	6	6	8	1	2	4	2	1	7	5	3	1	8	5	2	8	6	4	1	7	3	8	5	2	8	5	2	8	5	3	1	6	4	3	1	7	5	3	1	7	5	3	1	7	5	3 1	1	7	4	4	2	8	1	4	4	7	2	6	1	7	2	4	8	3	5	8	3	7	3	7	1	5	8	2	5	8	1	4	6	1	4	8	4	8	4	8	3	8	6	2	1	4	7	8	2	2	7	7	7	8	7	7	1	1	8	4	4	7	4	1	7	4	8	4	8	4	7	2	5	1	2	5	1	5	1	5	1	4	6	2	5	1	4	7	4	2	8	4	8	6	1	8	6	3	1	8	3	8	3	7];
w(w==1) = 0;
w(w==2) = 45;
w(w==3) = 90;
w(w==4) = 135;
w(w==5) = 180;
w(w==6) = 225;
w(w==7) = 270;
w(w==8) = 315;

% Control
y=[2 6 1 7 6 1 4 8 4 1 6 1 3 4 2 8	5	5	1	4	2	5	1	4	7	1	6	1	4	2	4	1	6	1	4	1	4	1	6	4	6	1	1	4	4	2	5	1	7	4	8	4	1	8	5	1	5	3	1	6	4	4	1	7	4	1	2	1	5	4	1	4	1	4	2	4	1	6	1	4	3	4	2	5	1	4	5	1	4	1	4	3	6	1	4 3	7	1	4	3	4	2	5	1	1	4	5	1	4	3	1	4	3	8	1	4	5	1	4	1	4	1	4	1	4	4	8	3	1	4	1	1	5	2	3	4	5	2	3	1	4	7	3	1	4	3	1	4	2	5	7	8	2	3	5	1	5	2	4	8	3	6	1	4	3	7	1	4	7	1	3	5	7	2	4];
y(y==1) = 0;
y(y==2) = 45;
y(y==3) = 90;
y(y==4) = 135;
y(y==5) = 180;
y(y==6) = 225;
y(y==7) = 270;
y(y==8) = 315;

% code for transposing a matrix
% averagematrix=[data(:,1:4)];
% average=reshape(averagematrix,[],1);

% Now convert degrees into radians 
alpha_deg=w';
alpha_rad = circ_ang2rad(alpha_deg); 

beta_deg =y';
beta_rad = circ_ang2rad(beta_deg); 

% Plot the circular directional heading 
figure(1)
subplot(2,2,1)
circ_plot(alpha_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Muscimol')
subplot(2,2,3)
circ_plot(alpha_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Muscimol')
subplot(2,2,2)
circ_plot(beta_rad,'pretty','bo',true,'linewidth',1,'color','r'),
title('Control')
subplot(2,2,4)
circ_plot(beta_rad,'hist',[],20,true,true,'linewidth',1,'color','r')
title('Control')

% code for the descriptive stats 
fprintf('\t\t\t\t\tALPHA\tBETA\n')

alpha_bar = circ_mean(alpha_rad);
beta_bar = circ_mean(beta_rad);

fprintf('Mean resultant vector:\t\t\t%.2f\t%.2f\n', circ_rad2ang([alpha_bar beta_bar]))

alpha_hat = circ_median(alpha_rad);
beta_hat = circ_median(beta_rad);

fprintf('Median:\t\t\t\t\t%.2f\t%.2f\n', circ_rad2ang([alpha_hat beta_hat]))

R_alpha = circ_r(alpha_rad);
R_beta = circ_r(beta_rad);

fprintf('R Length:\t\t\t\t%.2f\t%.2f\n',[R_alpha R_beta])

S_alpha = circ_var(alpha_rad);
S_beta = circ_var(beta_rad);

fprintf('Variance:\t\t\t\t%.2f\t%.2f\n',[S_alpha S_beta])

[s_alpha s0_alpha] = circ_std(alpha_rad);
[s_beta s0_beta] = circ_std(beta_rad);

fprintf('Standard deviation:\t\t\t%.2f\t%.2f\n',[s_alpha s_beta])
fprintf('Standard deviation 0:\t\t\t%.2f\t%.2f\n',[s0_alpha s0_beta])

b_alpha = circ_skewness(alpha_rad);
b_beta = circ_skewness(beta_rad);

fprintf('Skewness:\t\t\t\t%.2f\t%.2f\n',[b_alpha b_beta])

k_alpha = circ_kurtosis(alpha_rad);
k_beta = circ_kurtosis(beta_rad);

fprintf('Kurtosis:\t\t\t\t%.2f \t%.2f\n',[k_alpha k_beta])

fprintf('\n\n')

% code for the inferential stats 
fprintf('Inferential Statistics\n\nTests for Uniformity\n')
fprintf('\t\t\t\t\tALPHA\tBETA\n')

% Rayleigh test
p_alpha = circ_rtest(alpha_rad);
p_beta = circ_rtest(beta_rad);

fprintf('Rayleigh Test:\t\t\t\t%.2f\t%.2f\n',[p_alpha p_beta])
             
% Omnibus test
o_alpha = circ_otest(alpha_rad);
o_beta = circ_otest(beta_rad);

fprintf('Omnibus Test:\t\t\t\t %.2f\t%.2f\n',[o_alpha o_beta])

% Rao's spacing test
p_alpha = circ_raotest(alpha_rad);
p_beta = circ_raotest(beta_rad);

fprintf('Rao Spacing Test:\t\t\t%.2f\t%.2f\n',[p_alpha p_beta])

% V test
p_alpha = circ_vtest(alpha_rad,circ_ang2rad(0));
p_beta = circ_vtest(beta_rad,circ_ang2rad(0));

fprintf('V Test (r = 0):\t\t\t\t %.2f\t%.2f\n',[p_alpha p_beta])

% Watson-Williams multi-sample test for equal means
%   H0: the s populations have equal means
%   HA: the s populations have unequal means

% fprintf('\nWatson-Williams Multi-Sample test for equal means\n')
% [pval table] = circ_wwtest(alpha_rad,beta_rad)
% if (pval>.05)
%     fprintf('[\bpopulations have equal means]\b\n');
%     if (pval<.05)
%         fprintf('[\bpopulations DO NOT have equal means]\b\n');
%     end
% end
% fprintf('\nTests concerning Mean and Median angle\n')
% fprintf('\nDo not trust Cofidence means***********\n') % Had to create custom code for this
% fprintf('\nThey are probably wrong\n')
% % 95 percent confidence intervals for mean direction
% w = ones(size(alpha_rad));
% r = circ_r(alpha_rad,w,0,1);
% xi = 0.05;
% c2 = chi2inv((1-xi),1);
% n = sum(w,1);
% R = n.*r;
% t = zeros(size(r));
% i = 1:numel(r);
% sqrt(c2/2/n(i));
% t_alpha = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));
% 
% w = ones(size(beta_rad));
% r = circ_r(beta_rad,w,0,1);
% xi = 0.05;
% c2 = chi2inv((1-xi),1);
% n = sum(w,1);
% R = n.*r;
% t = zeros(size(r));
% i = 1:numel(r);
% sqrt(c2/2/n(i));
% t_beta = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));
% 
% %t_alpha = circ_confmean(alpha_rad,0.05);
% %t_beta = circ_confmean(beta_rad,0.05);
% 
% fprintf('Mean, up 95 perc. CI:\t\t\t %.2f\t%.2f\n', circ_rad2ang([alpha_bar+t_alpha beta_bar+t_beta]))
% fprintf('Mean, low 95 perc. CI:\t\t\t %.2f\t%.2f\n', circ_rad2ang([2*pi+alpha_bar-t_alpha 2*pi+beta_bar-t_beta])) 

