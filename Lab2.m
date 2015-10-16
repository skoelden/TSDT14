close all
clear all
constants

RyLoworder = R0*abs((1-a+a*exp(j*2*pi*theta))).^2;
ryLoworder = (1-a)*R0*not(tau)+a*R0*not(tau-1);

RyHighorder = R0*(rectpuls(theta/(2*theta0)) + rectpuls((1-theta)/(2*theta0)));
ryHighorder = 2*theta0*sinc(2*theta0*tau);

rySquare = ryHighorder(1)^2 + 2*ryHighorder.^2;
RySquare = ryHighorder(1)^2*not(theta) + 2*2*R0^2*theta0*(tripuls(theta/(4*theta0)) + tripuls((1-theta)/(4*theta0)));

ryRectified = ryHighorder(1)/(2*pi) + ryHighorder/4 + ryHighorder.^2/(4*pi*ryHighorder(1));
RyRectified = ryHighorder(1)./(2*pi)*not(theta) + R0/4*(rectpuls(theta/(2*theta0)) + rectpuls((1-theta)/(2*theta0))) + ...
    R0^2*2*theta0/(4*pi*ryHighorder(1))*(tripuls(theta/(4*theta0)) + tripuls((1-theta)/(4*theta0)));

RyAMSM = 1/4*(rectpuls((theta-omega0)/(2*theta0)) + rectpuls((theta+omega0)/(2*theta0))) + ...
    1/4*(rectpuls((1-theta+omega0)/(2*theta0)) + rectpuls((1-theta-omega0)/(2*theta0)));
    
squarednoise = idealfilterednoise.^2;

rectifiednoise = idealfilterednoise;
rectifiednoise(idealfilterednoise < 0) = 0;

AMSMnoise = idealfilterednoise.*cos(2*pi*omega0*n)';

[squaredACF, squaredPSD] = ACFe(squarednoise, 'bar', 's');
[rectifiedACF, rectifiedPSD] = ACFe(rectifiednoise, 'bar', 's');
[AMSMACF, AMSMPSD] = ACFe(AMSMnoise, 'bar', 's');

zeropoint = ceil(length(squaredACF)/2);

figure

subplot(1,3,1)
plot(0:1/(length(squaredPSD)-1):1, squaredPSD)
hold on
plot(0:1/(length(RySquare)-1):1, RySquare, 'r')
hold off
ylim([0, 0.4])

subplot(1,3,2)
plot(0:1/(length(rectifiedPSD)-1):1, rectifiedPSD)
hold on
plot(0:1/(length(RyRectified)-1):1, RyRectified, 'r')
hold off
ylim([0, 0.4])

subplot(1,3,3)
plot(0:1/(length(AMSMPSD)-1):1, AMSMPSD)
hold on
plot(0:1/(length(RyAMSM)-1):1, RyAMSM, 'r')
hold off
ylim([0, 0.4])