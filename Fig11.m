clear all
close all
clc
format longeng
depth=[30 35 40 45 50 55 60];
error_tangential=[2.1031 2.2985 2.6052 2.9974 3.5018 4.1026 4.8053];
error_radial=[1.8021 1.9963 2.2044 2.5026 2.8987 3.4034 3.9015];
plot(depth,error_tangential,'b-o')
hold on
plot(depth,error_radial,'r-d')
xlabel('Dipole depth (mm)')
ylabel('Localization error (mm)')
legend('Tangential','Radial')
title('Effect od dipole orientation on localization error')
grid on
axis([28 62 1.5 5])