close all
clear all
constants

RyHighorder = R0*(rectpuls(theta/(2*theta0)) + rectpuls((1-theta)/(2*theta0)));
ryHighorder = 2*theta0*sinc(2*theta0*tau);

alternatingnoise = idealfilterednoise;
for i = [1:length(alternatingnoise)]
    alternatingnoise(i) = alternatingnoise(i) * (-1)^(i+1);
end

decimatednoise = idealfilterednoise;
for i = [2:2:length(decimatednoise)-1]
    decimatednoise(i) = 0;
end

RyAlternating = rectpuls((theta-1/2)/(2*theta0));

RyDecimated = 1/4*(rectpuls(theta/(2*theta0)) + rectpuls((theta-1/2)/(2*theta0)) + rectpuls((1-theta)/(2*theta0)));

[alternatingACF, alternatingPSD] = ACFe(alternatingnoise, 'bar', 's');
[decimatedACF, decimatedPSD] = ACFe(decimatednoise, 'bar', 's');

figure

subplot(1,2,1)
plot(0:1/(length(alternatingPSD)-1):1, alternatingPSD)
hold on
plot(0:1/(length(RyAlternating)-1):1, RyAlternating, 'r')
hold off
%ylim([0, 0.4])

subplot(1,2,2)
plot(0:1/(length(decimatedPSD)-1):1, decimatedPSD)
hold on
plot(0:1/(length(RyDecimated)-1):1, RyDecimated, 'r')
hold off
%ylim([0, 0.4])
