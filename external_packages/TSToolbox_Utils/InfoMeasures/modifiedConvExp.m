function y = modifiedGaussConv(x,a)
Y = normpdf(-100:100,a(2),a(3));
x = a(1)*conv(x,Y,'same');
y=exp(x);