% test_MatdynModesToDebyeWallerFactor.m
clc;
clear;
close all;
%% main:
pwscf = ReadPwscfInput('tests/qe/si_cubic.1_scf.in');
mass = zeros(1, pwscf.system.nat);
for iType = 1 : pwscf.system.ntyp
    mass(pwscf.atomic_positions.types == pwscf.atomic_species.types(iType)) =...
        pwscf.atomic_species.masses(iType);
end

[qs, bands, eigenVecs] = ReadMatdynModes('tests/qe/matdyn.modes');

dwfs = MatdynModesToDebyeWallerFactor(qs, bands, eigenVecs, mass, 293);