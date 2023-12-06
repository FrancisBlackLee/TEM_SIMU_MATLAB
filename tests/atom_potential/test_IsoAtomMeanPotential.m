% test_IsoAtomMeanPotential.m
clc;
clear;
close all;
%% carbon film
carbonFilmDensity = 2.0e-24; % g / angs^3, see H.Lichte, 1998
carbonMolarMass = 12;
avogadroConst = 6.0221409e+23;
volPerAtom = carbonMolarMass / carbonFilmDensity / avogadroConst;
v0 = IsoAtomMeanPotential(volPerAtom, 6, 1);