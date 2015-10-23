close all
clear all
clc
constants

tic
RyLoworder = R0*abs((1-a+a*exp(-j*2*pi*theta))).^2;
ryLoworder = (a^2+(1-a)^2)*R0*not(tau)+a*(1-a)*R0*not(tau-1);

RyHighorder = R0*(rectpuls(theta/(2*theta0)) + rectpuls((1-theta)/(2*theta0)));
ryHighorder = 2*theta0*sinc(2*theta0*tau);

LOnoise = size(noise);
LOnoise(1) = (1-a)*noise(1);
for i = [2:length(noise)]
    LOnoise(i) = (1-a)*noise(i) + a*noise(i-1);
end

[bbutter, abutter] = butter(10, 2*theta0);
HOnoise = filter(bbutter, abutter, noise);

[loACF, loPSD] = ACFe(LOnoise, 'bar');
[hoACF, hoPSD] = ACFe(HOnoise, 'bar');

[loACFsm, loPSDsm] = ACFe(LOnoise, 'bar', 's');
[hoACFsm, hoPSDsm] = ACFe(HOnoise, 'bar', 's');
[hoACFoversm, hoPSDoversm] = ACFe(HOnoise, 'bar', 's', 11);
[hoACFundersm, hoPSDundersm] = ACFe(HOnoise, 'bar', 's', 3001 );

%Averaged periodogram
bins = 256;
len = length(LOnoise)/bins;
tmp = zeros([2*len-1 1])/bins;
for i = [0:bins-1]
    [tmpACF tmpPSD] = ACFe(LOnoise(i*len+1:(i+1)*len), 'bar');
    tmp = tmp + tmpPSD; 
end
loPSDav = tmp/bins;

bins = 256;
len = length(HOnoise)/bins;
tmp = zeros([2*len-1 1])/bins;
for i = [0:bins-1]
    [tmpACF tmpPSD] = ACFe(HOnoise(i*len+1:(i+1)*len), 'bar');
    tmp = tmp + tmpPSD; 
end
hoPSDav = tmp/bins;

bins = 2^12;
len = length(HOnoise)/bins;
tmp = zeros([2*len-1 1])/bins;
for i = [0:bins-1]
    [tmpACF tmpPSD] = ACFe(HOnoise(i*len+1:(i+1)*len), 'bar');
    tmp = tmp + tmpPSD; 
end
hoPSDoverav = tmp/bins;

bins = 16;
len = length(HOnoise)/bins;
tmp = zeros([2*len-1 1])/bins;
for i = [0:bins-1]
    [tmpACF tmpPSD] = ACFe(HOnoise(i*len+1:(i+1)*len), 'bar');
    tmp = tmp + tmpPSD; 
end
hoPSDunderav = tmp/bins;


toc
%%
fontSize = 16;

figure(1)
zeropoint = ceil(length(loACF)/2);
taus = -(length(loACF)-zeropoint):(length(loACF)-zeropoint);
plot(taus,loACF)
title('Bartletts estimate for low order filtered noise')
xlabel('Delay (samples)')
set(gca,'FontSize',fontSize)

figure(2)
zeropoint = ceil(length(loACF)/2);
stem([-20:20],loACF(zeropoint-20:zeropoint+20))
hold on
stem([-20:20],[fliplr(ryLoworder(2:length(tau))), ryLoworder(1:length(tau))],'rx'); 
hold off
title('Bartletts estimate for low order filtered noise')
xlabel('Delay (samples)')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(3)
zeropoint = ceil(length(hoACF)/2);
taus = -(length(hoACF)-zeropoint):(length(hoACF)-zeropoint);
plot(taus,loACF)
title('Bartletts estimate for high order filtered noise')
xlabel('Delay (samples)')
set(gca,'FontSize',fontSize)

figure(4)
zeropoint = ceil(length(hoACF)/2);
stem([-20:20],hoACF(zeropoint-20:zeropoint+20))
hold on
stem([-20:20],[fliplr(ryHighorder(2:length(tau))), ryHighorder(1:length(tau))],'rx'); 
hold off
title('Bartletts estimate for high order filtered noise')
xlabel('Delay (samples)')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(5)
plot(0:1/(length(loPSD)-1):1, loPSD)
hold on
plot(theta, RyLoworder, 'r')
hold off
title('PSD of low order filtered noise, raw')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(6)
plot(0:1/(length(loPSDsm)-1):1, loPSDsm)
hold on
plot(theta, RyLoworder, 'r')
hold off
title('PSD of low order filtered noise, smoothed')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(7)
plot(0:1/(length(loPSDav)-1):1, loPSDav)
hold on
plot(theta, RyLoworder, 'r')
hold off
title('PSD of low order filtered noise, averaged')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(8)
plot(0:1/(length(hoPSD)-1):1, hoPSD)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('PSD of high order filtered noise, raw')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(9)
plot(0:1/(length(hoPSDsm)-1):1, hoPSDsm)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('PSD of high order filtered noise, smoothed')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(10)
plot(0:1/(length(hoPSDav)-1):1, hoPSDav)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('PSD of high order filtered noise, averaged')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(11)
plot(0:1/(length(hoPSDoverav)-1):1, hoPSDoverav)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('Too much averageing')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(12)
plot(0:1/(length(hoPSDunderav)-1):1, hoPSDunderav)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('Not enough averageing')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(13)
plot(0:1/(length(hoPSDoversm)-1):1, hoPSDoversm)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('Too much smoothing')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)

figure(14)
plot(0:1/(length(hoPSDundersm)-1):1, hoPSDundersm)
hold on
plot(theta, RyHighorder, 'r')
hold off
title('Not enough smoothing')
xlabel('Normalized frequency, \theta')
legend('Estimated', 'Theoretical')
set(gca,'FontSize',fontSize)