close all
clear all
constants
fontSize = 10.5;

RyHighorder = R0*(rectpuls(theta/(2*theta0)) + rectpuls((1-theta)/(2*theta0)));
ryHighorder = 2*theta0*sinc(2*theta0*tau);

%Alternating noise
alternatingnoise = idealfilterednoise;
for i = [1:length(alternatingnoise)]
    alternatingnoise(i) = alternatingnoise(i) * (-1)^(i+1);
end

%Decimated noise
decimatednoise = idealfilterednoise;
for i = [2:2:length(decimatednoise)-1]
    decimatednoise(i) = 0;
end

%Fourier transform of the two noises
RyAlternating = rectpuls((theta-1/2)/(2*theta0));

RyDecimated = 1/4*(rectpuls(theta/(2*theta0)) + rectpuls((theta-1/2)/(2*theta0)) + rectpuls((1-theta)/(2*theta0)));

% Calculating the PSD of the signals
[alternatingACF, alternatingPSD] = ACFe(alternatingnoise, 'bar', 's');
[decimatedACF, decimatedPSD] = ACFe(decimatednoise, 'bar', 's');

% Calculating the PSD of the signals
[rawalternatingACF, rawalternatingPSD] = ACFe(alternatingnoise, 'bar');
[rawdecimatedACF, rawdecimatedPSD] = ACFe(decimatednoise, 'bar');


%% plotting

figure(1)
plot(0:1/(length(alternatingPSD)-1):1, alternatingPSD)
hold on
plot(0:1/(length(RyAlternating)-1):1, RyAlternating, 'r')
hold off
legend('Estimate', 'Theoretical'), title('PSD of alternating signal')
xlabel('Normalized frequency, \theta')
set(gca,'FontSize',fontSize)

figure(2)
plot(0:1/(length(decimatedPSD)-1):1, decimatedPSD)
hold on
plot(0:1/(length(RyDecimated)-1):1, RyDecimated, 'r')
hold off
legend('Estimate', 'Theoretical'), title('PSD of decimated signal')
xlabel('Normalized frequency, \theta')
set(gca,'FontSize',fontSize)

