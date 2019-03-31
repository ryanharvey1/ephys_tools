function [I1,I2]=linecirc_intersect(P1,P2,C,r)
  A = P1-C;
  B = P2-P1;
  d2 = dot(B,B);
  t = (r^2-dot(A,A))*d2+dot(A,B)^2;
  Q = P1-dot(A,B)/d2*B;
  t2 = sqrt(t)/d2*B;
  I1 = Q + t2;
  I2 = Q - t2;
end