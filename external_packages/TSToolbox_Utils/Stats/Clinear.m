function [R rdw rup] = Clinear(linVar,theta)

nbIterations = 1000;

n = length(theta);
C = cos(theta);
S = sin(theta);

r12 = corrcoef([linVar C]);
r12 = r12(1,2);
r13 = corrcoef([linVar S]);
r13 = r13(1,2);
r23 = corrcoef([C S]);
r23 = r23(1,2);


R = (r12^2+r13^2-2*r12*r13*r23)/(1-r23^2);

Rrand = [];

for i=1:nbIterations

	t = theta(randperm(n));
	l = linVar(randperm(n));
	
	C = cos(t);
	S = sin(t);

	r12 = corrcoef([linVar C]);
	r12 = r12(1,2);
	r13 = corrcoef([linVar S]);
	r13 = r13(1,2);
	r23 = corrcoef([C S]);
	r23 = r23(1,2);
	
	r = (r12^2+r13^2-2*r12*r13*r23)/(1-r23^2);
	Rrand = [Rrand, r];


end


Rrand = sort(Rrand);
rdw = Rrand(round(0.05*nbIterations));
rup = Rrand(round(0.95*nbIterations));
