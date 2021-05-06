% test_ElectronChargeDensity.m
clc;
clear;
close all;

%% main
r = linspace(0.001, 1.0, 1000);
rho = ElectronDensity(32, r);
rho = rho * 4 * pi .* r.^2;

plot(r, rho);
% data = [r;rho]';
% save('tests\charge_density_Ge.txt', 'data', '-ascii', '-double', '-tabs');