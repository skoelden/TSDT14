close all
clear all
constants

RyLoworder = R0*abs((1-a+a*exp(j*2*pi*theta))).^2;
ryLoworder = (1-2*a+a^2)*R0*not(tau)+(1-a)*a*R0*not(tau-1);

RyHighorder = R0*(rectpuls(theta/(2*theta0)) + rectpuls((1-theta)/(2*theta0)));
ryHighorder = 2*theta0*sinc(2*theta0*tau);

filterednoise = size(noise);
filterednoise(1) = (1-a)*noise(1);
for i = [2:length(noise)]
    filterednoise(i) = (1-a)*noise(i) + a*noise(i-1);
end

[bbutter, abutter] = butter(30, 2*theta0);
idealfilterednoise = filter(bbutter, abutter, noise);

[filteredACF, filteredPSD] = ACFe(filterednoise, 'bar', 's');
[idealfilteredACF, idealfilteredPSD] = ACFe(idealfilterednoise, 'bar', 's');

%Averaged periodogram
bins = 16;
len = length(idealfilterednoise)/bins;
tmp = zeros([2*len-1 1])/bins;
for i = [0:bins-1]
    tmp = tmp + ACFe(idealfilterednoise(i*len+1:(i+1)*len), 'bar');
end
idealfilteredACFavg = tmp/bins;
idealfilteredPSDavg = abs(fft(idealfilteredACFavg));

figure(3)
subplot(2,2,1)
plot(0:1/(length(filteredPSD)-1):1, real(filteredPSD))
hold on
plot(theta, RyLoworder, 'r')
hold off
subplot(2,2,2)
plot(0:1/(length(idealfilteredPSD)-1):1, real(idealfilteredPSD))
hold on
plot(theta, RyHighorder, 'r')
hold off
subplot(2,2,3)
zeropoint = ceil(length(idealfilteredACF)/2);
stem(tau,filteredACF(zeropoint:zeropoint+(length(tau)-1)));
hold on
stem(tau, ryLoworder, 'r')
hold off
subplot(2,2,4)
zeropoint = ceil(length(idealfilteredACF)/2);
plot(tau,idealfilteredACF(zeropoint:zeropoint+(length(tau)-1)));
hold on
plot(tau, ryHighorder, 'r')
hold off

%plot PSD with low order filtering
figure(1),
plot(0:1/(length(filteredPSD)-1):1, real(filteredPSD) , theta, RyLoworder , 'r' )
legend('Estimated PSD','Theoretical PSD')
title('Comparison of estimated PSD and theoretical PSD filtered with low order filter')
