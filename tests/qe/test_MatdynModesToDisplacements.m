% test_MatdynModesToDisplacements.m
clc;
clear;
close all;
%% parse pwscf and matdyn.modes
pwscf = ReadPwscfInput('tests/qe/si_cubic.1_scf.in');
[unitCell, a1, a2, a3, b1, b2, b3] = PwscfInputToFracCoords(pwscf);
mass = zeros(1, pwscf.system.nat);
for iType = 1 : pwscf.system.ntyp
    mass(pwscf.atomic_positions.types == pwscf.atomic_species.types(iType)) =...
        pwscf.atomic_species.masses(iType);
end

[qs, bands, eigenVecs] = ReadMatdynModes('tests/qe/matdyn.modes');

%% calculate displacements:
nConfig = 2;
superCell = TileUnitCell(unitCell, [8, 8, 8]);
convMat = BasesToConvMat(a1, a2, a3);

displacements = MatdynModesToDisplacements(superCell, a1, a2, a3, b1, b2, b3,...
    qs, bands, eigenVecs, mass, 293, nConfig);

coords = convMat * superCell(3 : 5, :);

for iConfig = 1 : nConfig
    tmpCoords = coords + displacements(:, :, iConfig);
    figure;
    scatter3(tmpCoords(1, :), tmpCoords(2, :), tmpCoords(3, :), '.');
    title(['config. ', num2str(iConfig)]);
end