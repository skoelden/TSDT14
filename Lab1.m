close all
clear all
clc
constants
tic
RyLoworder = R0*abs((1-a+a*exp(-j*2*pi*theta))).^2;
ryLoworder = (a^2+(1-a)^2)*R0*not(tau)+a*(1-a)*R0*not(tau-1);

RyHighorder = R0*(rectpuls(theta/(2*theta0)) + rectpuls((1-theta)/(2*theta0)));
ryHighorder = 2*theta0*sinc(2*theta0*tau);

filterednoise = size(noise);
filterednoise(1) = (1-a)*noise(1);
for i = [2:length(noise)]
    filterednoise(i) = (1-a)*noise(i) + a*noise(i-1);
end

[bbutter, abutter] = butter(10, 2*theta0);
idealfilterednoise = filter(bbutter, abutter, noise);

[BARfilteredACF, BARfilteredPSD] = ACFe(filterednoise, 'bar');
[BARidealfilteredACF, BARidealfilteredPSD] = ACFe(idealfilterednoise, 'bar');

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
fontSize = 16;

figure(1)
zeropoint = ceil(length(BARfilteredACF)/2);
taus = -(length(BARfilteredACF)-zeropoint):(length(BARfilteredACF)-zeropoint);
plot(taus,BARfilteredACF)
title('Bartletts estimate for low order filtered noise')
xlabel('Delay (samples)')
set(gca,'FontSize',fontSize)

figure(2)
zeropoint = ceil(length(BARfilteredACF)/2);
stem([-20:20],BARfilteredACF(zeropoint-20:zeropoint+20))
hold on
stem([-20:20],[fliplr(ryLoworder(2:length(tau))), ryLoworder(1:length(tau))],'rx'); 
hold off
title('Bartletts estimate for low order filtered noise')
xlabel('Delay (samples)')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(3)
zeropoint = ceil(length(BARidealfilteredACF)/2);
taus = -(length(BARidealfilteredACF)-zeropoint):(length(BARidealfilteredACF)-zeropoint);
plot(taus,BARfilteredACF)
title('Bartletts estimate for high order filtered noise')
xlabel('Delay (samples)')
set(gca,'FontSize',fontSize)

figure(4)
zeropoint = ceil(length(BARidealfilteredACF)/2);
stem([-20:20],BARidealfilteredACF(zeropoint-20:zeropoint+20))
hold on
stem([-20:20],[fliplr(ryHighorder(2:length(tau))), ryHighorder(1:length(tau))],'rx'); 
hold off
title('Bartletts estimate for high order filtered noise')
xlabel('Delay (samples)')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(5)
plot(0:1/(length(BARfilteredPSD)-1):1, BARfilteredPSD)
hold on
plot(theta, RyLoworder, 'r')
hold off
title('PSD of low order filtered noise, raw')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(6)
plot(0:1/(length(BARsmoothedfilteredPSD)-1):1, BARsmoothedfilteredPSD)
hold on
plot(theta, RyLoworder, 'r')
hold off
title('PSD of low order filtered noise, smoothed')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(7)
plot(0:1/(length(BARaveragedfilteredPSD)-1):1, BARaveragedfilteredPSD)
hold on
plot(theta, RyLoworder, 'r')
hold off
title('PSD of low order filtered noise, averaged')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(8)
plot(0:1/(length(BARidealfilteredPSD)-1):1, BARidealfilteredPSD)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('PSD of high order filtered noise, raw')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(9)
plot(0:1/(length(BARsmoothedidealfilteredPSD)-1):1, BARsmoothedidealfilteredPSD)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('PSD of high order filtered noise, smoothed')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(10)
plot(0:1/(length(BARaveragedidealfilteredPSD)-1):1, BARaveragedidealfilteredPSD)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('PSD of high order filtered noise, averaged')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)