% test_RadialAtomPotential.m
clc;
clear;
close all;
%% potential:
r = linspace(0.001, 1.0, 1000);
radPot = RadialAtomPotential(32, r);

plot(r, radPot);

data = [r; radPot]';
save('tests\radial_potential_Ge.txt', 'data', '-ascii', '-double', '-tabs');