function [rets] = gennetvar(response, SNR, numvar)

if nargin < 1
    response = 0.4;
    SNR = 10;
    numvar = 3;
end

a = -1;
b = 1;
coeffvec = a + (b-a).*rand(numvar,1);
errorterm = response/SNR*randn(1,1);
rets = zeros(numvar,1);

rets = (response - errorterm)/sum(coeffvec) * coeffvec;
