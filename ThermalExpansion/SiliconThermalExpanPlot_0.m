% SiliconThermalExpanPlot_0.m
% Ref: https://trc.nist.gov/cryogenics/materials/Silicon/Silicon.htm
clc;
close all;
clear all;
%% main:
T = 0 : 1 : 600;

a = 1.00500E-05;
b = -5.99688E-06;
c = 1.25574E-06;
d = -1.12086E-07;
e = 3.63225E-09;
f = 2.67708E-02;
g = -1.22829E-04;
h = 1.62544E-18;
i = 4.72374E+02;
j = -3.58796E+04;
k = -1.24191E+07;
l = 1.25972E+09;

expanCoeff = (4.8e-5 * T.^3 + (a * T.^5 + b * T.^5.5 + c * T.^6 ...
    + d * T.^6.5 + e * T.^7) .* ((1 + erf(T - 15)) / 2)) ...
    .* ((1 - erf(0.2 * (T - 52))) / 2) + ((-47.6 + f * (T - 76).^2 ...
    + g * (T - 76).^3 + h * (T - 76).^9) .* ((1 + erf(0.2 * (T - 52))) / 2)) ...
    .* ((1 - erf(0.1 * (T - 200))) / 2) + ((i + j ./ T + k ./ T.^2 ...
    + l ./ T.^3) .* ((1 + erf(0.1 * (T - 200))) / 2));

figure;
plot(T, expanCoeff);