function y = modifiedGaussConv(x,a)

Y = normpdf(-100:100,0,a(1).^2);
x = conv(x,Y,'same');
y = exp(a(2)*x);