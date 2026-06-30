clear all
close all
clc
format longeng
depth=[30 35 40 45 50 55 60];
Three_shell=[4.2 4.8273 5.6136 6.5 7.5091 8.5750 9.8];
proposed=[1.5 1.7136 2 2.4136 2.9273 3.3409 3.8];
plot(depth,proposed,'b-o')
hold on
plot(depth,Three_shell,'r--s')
xlabel('Dipole depth (mm)')
ylabel('Amplitude bias (%)')
title('Amplitude bias versus dipole depth')
grid on
axis([28 62 1 11])
legend('Propose model','3-shell model')