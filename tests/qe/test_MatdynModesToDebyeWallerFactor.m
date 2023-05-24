% test_MatdynModesToDebyeWallerFactor.m
clc;
clear;
close all;
%% main:
pwscf = ReadPwscfInput('tests/qe/si_cubic.1_scf.in');
pwscf = ExpandPwscfAtomMass(pwscf);

[qs, bands, eigenVecs] = ReadMatdynModes('tests/qe/matdyn.modes');

dwfs = MatdynModesToDebyeWallerFactor(qs, bands, eigenVecs, pwscf.mass, 293);