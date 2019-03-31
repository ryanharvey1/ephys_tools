function Y = KalmanFilter(X,hgamma)
%Input
%X:         time series data
%hgamma:    hyperparameters
%Output
%Y:         estimates
%
%version 0.1
%2010-3-31 Takeaki Shimokawa
%Copyright (c) 2010, Takeaki Shimokawa All rights reserved.

T = diff(X);    %ISIs
N = length(T);  %number of ISIs

mu = mean(T);
Lv = 0;         %local variation
for i = 1 : N-1
    Lv = Lv + ((T(i)-T(i+1))./(T(i)+T(i+1))).^2;
end
Lv = Lv.*3./(N-1);

EL = zeros(2,N);        %[lambda(i|i-1), lambda(i|i)]
VL = zeros(2,N);        %variance of lambda
EL_N = zeros(1,N);      %lambda(i|N)
VL_N = zeros(1,N);
COVL_N = zeros(1,N);    %V(1,1)(i+1,i|N)

EK = zeros(2,N);        %[kappa(i!i-1), kappa(i|i)]
VK = zeros(2,N);        %variance of kappa
EK_N = zeros(1,N);      %kappa(i|N)
VK_N = zeros(1,N);
COVK_N = zeros(1,N);    %V(2,2)(i+1,i|N)

VLK = zeros(2,N);       %covariance of lambda and kappa
VLK_N = zeros(1,N);

for i = 1 : N
    %---prediction
    if i==1
        EL(1,i) = 1 ./ mu;
        VL(1,i) = (EL(1,i)./2).^2;
        EK(1,i) = (3./Lv-1)./2;
        VK(1,i) = (EK(1,i)./2).^2;
        VLK(1,i) = 0;
    else
        EL(1,i) = EL(2,i-1);
        EK(1,i) = EK(2,i-1);
        VL(1,i) = VL(2,i-1) + T(i-1).*hgamma(1).^2;
        VK(1,i) = VK(2,i-1) + T(i-1).*hgamma(2).^2;
        VLK(1,i) = VLK(2,i-1);
    end
    %---filtering
    detV = VL(1,i).*VK(1,i) - VLK(1,i).^2;
    alpha = VK(1,i) ./ detV;
    beta = - VLK(1,i) ./ detV;
    gamma = VL(1,i) ./ detV;
    EL1 = EL(1,i);
    EK1 = EK(1,i);
    EL2 = EL1;
    EK2 = EK1;
    for m = 1 : 3
        %filtering lambda
        A = alpha.*EL1 - beta.*(EK2-EK1) - EK2.*T(i);
        EL2 = (A+sqrt(A.*A+4.*alpha.*EK2))./(2.*alpha);
        %filtering kappa --- using the bisection method
        mid = EK2;
        f = - gamma.*(mid-EK1) - beta.*(EL2-EL1) ...
            + (1+log(EL1.*T(i))-EL1.*T(i)) + log(mid) - psi(mid);
        if f >= 0
            right = mid;
            while f >= 0
                left = right;
                right = left.*2;
                f = - gamma.*(right-EK1) - beta.*(EL2-EL1) ...
                    + (1+log(EL1.*T(i))-EL1.*T(i)) + log(right) - psi(right);
            end
        elseif f < 0
            left = mid;
            while f < 0
                right = left;
                left = right./2;
                f = - gamma.*(left-EK1) - beta.*(EL2-EL1) ...
                    + (1+log(EL1.*T(i))-EL1.*T(i)) + log(left) - psi(left);
            end
        end
        for j = 1 : 10
            mid = (left + right) ./ 2;
            f = - gamma.*(mid-EK1) - beta.*(EL2-EL1) ...
                + (1+log(EL1.*T(i))-EL1.*T(i)) + log(mid) - psi(mid);
            if f >= 0
                left = mid;
            elseif f < 0
                right = mid;
            end
        end
        EK2 = (left + right) ./ 2;
    end
    EL(2,i) = EL2;
    EK(2,i) = EK2;
    a = alpha + EK2 ./ EL2.^2;
    b = beta - 1 ./ EL2 + T(i);
    c = gamma - 1 ./ EK2 + psi(1,EK2);
    detV = a.*c - b.^2;
    VL(2,i) = c ./ detV;
    VLK(2,i) = - b ./ detV;
    VK(2,i) = a ./ detV;
end
%---smoothing
EL_N(N) = EL(2,N);
VL_N(N) = VL(2,N);
EK_N(N) = EK(2,N);
VK_N(N) = VK(2,N);
VLK_N(N) = VLK(2,N);
for i = N-1 : -1 : 1
    Etheta1 = [EL(2,i), EK(2,i)]';
    Etheta2 = [EL(1,i+1), EK(1,1+i)]';
    EthetaN = [EL_N(i+1), EK_N(i+1)]';
    Vtheta1 = [VL(2,i), VLK(2,i); VLK(2,i), VK(2,i)];
    Vtheta2 = [VL(1,i+1), VLK(1,i+1); VLK(1,i+1), VK(1,i+1)];
    VthetaN = [VL_N(i+1), VLK_N(i+1); VLK_N(i+1), VK_N(i+1)];
    
    H = Vtheta1 / Vtheta2;
    Etheta = Etheta1 + H * ( EthetaN - Etheta2 );
    Vtheta = Vtheta1 + H * ( VthetaN - Vtheta2 ) * H';
    COVtheta = H * VthetaN;
    
    EL_N(i) = Etheta(1);
    EK_N(i) = Etheta(2);
    VL_N(i) = Vtheta(1,1);
    VK_N(i) = Vtheta(2,2);
    VLK_N(i) = Vtheta(1,2);
    COVL_N(i) = COVtheta(1,1);
    COVK_N(i) = COVtheta(2,2);
end
Y = [EL_N; VL_N; COVL_N; EK_N; VK_N; COVK_N];
