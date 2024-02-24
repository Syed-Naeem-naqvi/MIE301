clear;

% Practice 1
a = 5;
b = 1:10;
b = b';

clear all

c = 1:40;
clc;
last_ten = c(31:end);
time = 0:0.1:4;
size = length(time);
figure
hold on
plot(time, sin(time), 'ro')
title('Test Graph', FontSize=15)
xlabel('Time', FontSize=15)
ylabel('sin(time)', fontsize=15)
grid

plot(time, cos(time), 'x')
plot(time, tan(time),'*', Color='red')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2, 1, 1);
x = linspace(0, 10);
y1 = sin(x);
plot(x, y1)

subplot(2, 1, 2);
y2 = sin(5*x);
plot(x, y2)