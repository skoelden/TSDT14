close all
clear all
constants
tic
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

[BARfilteredACF, BARfilteredPSD] = ACFe(filterednoise, 'bar');

[BARsmoothedfilteredACF, BARsmoothedfilteredPSD] = ACFe(filterednoise, 'bar', 's');
[BARsmoothedidealfilteredACF, BARsmoothedidealfilteredPSD] = ACFe(idealfilterednoise, 'bar', 's');

%Averaged periodogram
bins = 256;
len = length(filterednoise)/bins;
tmp = zeros([2*len-1 1])/bins;
for i = [0:bins-1]
    [tmpACF tmpPSD] = ACFe(filterednoise(i*len+1:(i+1)*len), 'bar');
    tmp = tmp + tmpPSD; 
end
BARaveragedfilteredPSD = tmp/bins;

bins = 256;
len = length(idealfilterednoise)/bins;
tmp = zeros([2*len-1 1])/bins;
for i = [0:bins-1]
    [tmpACF tmpPSD] = ACFe(idealfilterednoise(i*len+1:(i+1)*len), 'bar');
    tmp = tmp + tmpPSD; 
end
BARaveragedidealfilteredPSD = tmp/bins;

toc
%%

figure(1)
zeropoint = ceil(length(BARfilteredACF)/2);
taus = -(length(BARfilteredACF)-zeropoint):(length(BARfilteredACF)-zeropoint);
plot(taus,BARfilteredACF)
title('Bartletts estimate of low order filtered noise')
xlabel('Delay (samples)')
print('TSDT14/Images/BARloACF', '-dpng')

figure(2)
zeropoint = ceil(length(BARfilteredACF)/2);
stem(tau,BARfilteredACF(zeropoint:zeropoint+20))
hold on
stem(tau,ryLoworder(1:length(tau)),'rx'); 
hold off
title('Bartletts estimate of low order filtered noise')
xlabel('Delay (samples)')
print('TSDT14/Images/BTloACFfirst20', '-dpng')

figure(3)
plot(0:1/(length(BARfilteredPSD)-1):1, BARfilteredPSD)
hold on
plot(theta, RyLoworder, 'r')
hold off
title('PSD from Bartletts estimate')
xlabel('Normalized frequency, \theta')
print('TSDT14/Images/BARloACF', '-dpng')

figure(4)
plot(0:1/(length(BARsmoothedfilteredPSD)-1):1, BARsmoothedfilteredPSD)
hold on
plot(theta, RyLoworder, 'r')
hold off
title('PSD from Bartletts estimate, smoothed')
xlabel('Normalized frequency, \theta')
print('TSDT14/Images/BTloPSDsmoothed', '-dpng')

figure(5)
plot(0:1/(length(BARaveragedfilteredPSD)-1):1, BARaveragedfilteredPSD)
hold on
plot(theta, RyLoworder, 'r')
hold off
title('PSD from Bartletts estimate, averaged')
xlabel('Normalized frequency, \theta')
print('TSDT14/Images/BARloPSDsmoothed', '-dpng')

figure(6)
plot(0:1/(length(BARsmoothedidealfilteredPSD)-1):1, BARsmoothedidealfilteredPSD)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('PSD from Bartletts estimate, smoothed')
xlabel('Normalized frequency, \theta')
print('TSDT14/Images/BARhoPSDsmoothed', '-dpng')

figure(7)
plot(0:1/(length(BARaveragedidealfilteredPSD)-1):1, BARaveragedidealfilteredPSD)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('PSD from Bartletts estimate, averaged')
xlabel('Normalized frequency, \theta')
print('TSDT14/Images/BARhoPSDsmoothed', '-dpng')
