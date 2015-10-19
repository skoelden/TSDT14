function [acf per] = ACFe(in, type, option)
L = 101;
N = length(in);
acf = zeros(N,1);

for kloop = 1:N
    k = kloop - 1;
    tmp = 0;
    for n = 1:(N-k)
        tmp = tmp + in(n)*in(n+k);
    end
    
    if strcmp(type, 'bt')
        acf(kloop) = 1/(N-k) * tmp;
    elseif strcmp(type, 'bar')
        acf(kloop) = 1/(N) * tmp;
    else
        error('No correct type specified')
    end
end

acf2 = zeros(2*length(acf)-1, 1);
len = length(acf);
for i = 1:len
    acf2(i) = acf(len+1-i);
    acf2(2*len-i) = acf(len+1-i);
end

acf = acf2;

if ~exist('option', 'var')
    per = abs(fft(acf));
    
elseif strcmp(option, 's')
    midway = ceil(size(acf2)/2);
    ham = zeros(size(acf2));
    ham(midway-floor(L/2):midway+floor(L/2)) = hamming(L);
    acf2 = acf2.* ham;
    per = abs(fft(acf2));
end
end