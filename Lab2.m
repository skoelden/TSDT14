close all
clear all
constants

%filter time discrete
RyLoworder = R0*abs((1-a+a*exp(1j*2*pi*theta))).^2;
ryLoworder = (1-a)*R0*not(tau)+a*R0*not(tau-1);

%Filter, fourier transformed
RyHighorder = R0*(rectpuls(theta/(2*theta0)) + rectpuls((1-theta)/(2*theta0)));
ryHighorder = 2*theta0*sinc(2*theta0*tau);

% Input and its fourier transform
rySquare = ryHighorder(1)^2 + 2*ryHighorder.^2;
RySquare = ryHighorder(1)^2*not(theta) + 2*2*R0^2*theta0*(tripuls(theta/(4*theta0)) + tripuls((1-theta)/(4*theta0)));

ryRectified = ryHighorder(1)/(2*pi) + ryHighorder/4 + ryHighorder.^2/(4*pi*ryHighorder(1));
RyRectified = ryHighorder(1)./(2*pi)*not(theta) + R0/4*(rectpuls(theta/(2*theta0)) + rectpuls((1-theta)/(2*theta0))) + ...
    R0^2*2*theta0/(4*pi*ryHighorder(1))*(tripuls(theta/(4*theta0)) + tripuls((1-theta)/(4*theta0)));

RyAMSM = 1/4*(rectpuls((theta-omega0)/(2*theta0)) + rectpuls((theta+omega0)/(2*theta0))) + ...
    1/4*(rectpuls((1-theta+omega0)/(2*theta0)) + rectpuls((1-theta-omega0)/(2*theta0)));

% filtered systemnoise
squarednoise = idealfilterednoise.^2;

rectifiednoise = idealfilterednoise;
rectifiednoise(idealfilterednoise < 0) = 0;

AMSMnoise = idealfilterednoise.*cos(2*pi*omega0*n)';

% ACF, PSD
[squaredACF, squaredPSD] = ACFe(squarednoise, 'bar', 's');
[rectifiedACF, rectifiedPSD] = ACFe(rectifiednoise, 'bar', 's');
[AMSMACF, AMSMPSD] = ACFe(AMSMnoise, 'bar', 's');

% raw estimates
[rawsquaredACF, rawsquaredPSD] = ACFe(squarednoise, 'bar');
[rawrectifiedACF, rawrectifiedPSD] = ACFe(rectifiednoise, 'bar');
[rawAMSMACF, rawAMSMPSD] = ACFe(AMSMnoise, 'bar');
%zeropoint = ceil(length(squaredACF)/2);

%% raw estimates of the PSDs
fontSize = 16;

figure (1)
plot(0:1/(length(rawsquaredPSD)-1):1, rawsquaredPSD)
hold on
plot(0:1/(length(RySquare)-1):1, RySquare, 'r')
hold off
ylim([0, 0.5])
legend('Estimate', 'Theoretical'), title ('PSD of squarer, raw vs. theoretical ')
xlabel('Normalized frequency, \theta')
set(gca,'FontSize',fontSize)

figure (2)
plot(0:1/(length(rawrectifiedPSD)-1):1, rawrectifiedPSD)
hold on
plot(0:1/(length(RyRectified)-1):1, RyRectified, 'r')
hold off
ylim([0, 0.4]),legend('Estimate', 'Theoretical'), title ('PSD of half-wave, raw vs. theoretical')
xlabel('Normalized frequency, \theta')
set(gca,'FontSize',fontSize)

figure (3)
plot(0:1/(length(rawAMSMPSD)-1):1, rawAMSMPSD)
hold on
plot(0:1/(length(RyAMSM)-1):1, RyAMSM, 'r')
hold off
ylim([0, 0.4]),legend('Estimate', 'Theoretical'), title ('PSD of Am-SC modulator, raw vs. theoretical')
xlabel('Normalized frequency, \theta')
set(gca,'FontSize',fontSize)

%% Smoothed estimates of the PSDs
figure (4)
plot(0:1/(length(squaredPSD)-1):1, squaredPSD)
hold on
plot(0:1/(length(RySquare)-1):1, RySquare, 'r')
hold off
ylim([0, 0.6]), legend('Estimate', 'Theoretical');
title('PSD of squarer, smoothed vs. theoretical')
xlabel('Normalized frequency, \theta')
set(gca,'FontSize',fontSize)

figure(5)
plot(0:1/(length(rectifiedPSD)-1):1, rectifiedPSD)
hold on
plot(0:1/(length(RyRectified)-1):1, RyRectified, 'r')
hold off
ylim([0, 0.4]), legend('Estimate', 'Theoretical');
title('PSD of half-wave, smoothed vs. theoretical')
xlabel('Normalized frequency, \theta')
set(gca,'FontSize',fontSize)

figure(6)
plot(0:1/(length(AMSMPSD)-1):1, AMSMPSD)
hold on
plot(0:1/(length(RyAMSM)-1):1, RyAMSM, 'r')
hold off
ylim([0, 0.4]), legend('Estimate', 'Theoretical');
title('PSD of AM-SC, smoothed vs. theoretical')
xlabel('Normalized frequency, \theta')
set(gca,'FontSize',fontSize)

figure(7)
plot(0:1/(length(RyHighorder)-1):1, RyHighorder)
ylim([0, 0.4]), 
title('PSD of the input')
xlabel('Normalized frequency, \theta')
set(gca,'FontSize',fontSize)
%% Histograms of the outputs
figure(8)
h1 = histogram(noise);
h1.Normalization = 'probability';
title('The histogram of the input')
xlabel('Vaule'), ylabel('Probability')
set(gca,'FontSize',fontSize)

figure(9)
h2 = histogram(squarednoise);
h2.Normalization = 'probability';
title('Histogram of the output, squarer')
xlabel('Vaule'), ylabel('Probability')
set(gca,'FontSize',fontSize)

figure(10)
h3=histogram(rectifiednoise);
h3.Normalization = 'probability';
title('Histogram of the output, halfwave rectifier')
xlabel('Vaule'), ylabel('Probability')
set(gca,'FontSize',fontSize)

figure(11)
h4 = histogram(AMSMnoise);
h4.Normalization = 'probability';
title('Histogram of the output, AM-SC')
xlabel('Vaule'), ylabel('Probability')
set(gca,'FontSize',fontSize)

