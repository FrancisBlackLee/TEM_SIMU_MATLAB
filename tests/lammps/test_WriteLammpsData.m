% test_WriteLammpsData.m
clc;
clear;
close all;
%% main:
[unitCell, convMat, cutoff] = CrystalAdvisor('tests/lammps/AlGaAs2_mp-1228891_symmetrized.cif', [1, 1, 0], true);

unitCell = RemoveSymmetricAtoms(unitCell);

figure;
PlotUnitCell3D(convMat, unitCell);

nc = [3, 1, 2];
supercell = TileUnitCell(unitCell, nc);

figure;
PlotUnitCell3D(convMat, supercell);

supercell(3 : 5, :) = convMat * supercell(3 : 5, :);

supercellConvMat = [nc(1) * convMat(:, 1), nc(2) * convMat(:, 2), nc(3) * convMat(:, 3)];

WriteLammpsData('tests/lammps/AlGaAs2_3x1x2.lammps', supercellConvMat, supercell([1, 3 : 5], :));