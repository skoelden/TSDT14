N = 2^16;
n = linspace(0,N,N);
Ts = 1;
fs = N/Ts;
T = Ts/N;

R0 = 1;
theta = [0:0.001:1];

a = 0.4;
%n = 0:10;

tau = 0:50;

theta0 = 0.1;
rect1 = [-0.5:0.00001:0.5]<theta0;
rect2 = [-0.5:0.00001:0.5]>-theta0;
rect = rect1 .* rect2;

noise = randn(N,1);
[bbutter, abutter] = butter(30, 2*theta0);
idealfilterednoise = filter(bbutter, abutter, noise);
omega0 = 0.2;