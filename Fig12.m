clear all
close all
clc
format longeng
depth=[30 35 40 45 50 55 60];
Three_shell=[3.51916 4.00697 4.80139 5.93031 7.14286 8.53659 10.0418];
proposed=[2.09756 2.30662 2.5993 3.00348 3.50523 4.10453 4.81533];
plot(depth,proposed,'b-o')
hold on
plot(depth,Three_shell,'r--s')
xlabel('Dipole depth (mm)')
ylabel('Localization error (mm)')
title('Localization error versus dipole depth (tangential)')
grid on
axis([28 62 1 11])
legend('Propose anisotropic-eccentric model','3-shell isotropic concentric model')