% test_WriteLammpsData.m
clc;
clear;
close all;
%% main:
[unitCell, convMat, cutoff] = CrystalAdvisor('tests/lammps/AlGaAs2_mp-1228891_symmetrized.cif', [1, 1, 0], true);

unitCell = RemoveSymmetricAtoms(unitCell);

figure;
PlotUnitCell3D(convMat, unitCell);

supercell = TileUnitCell(unitCell, [3, 1, 2]);

figure;
PlotUnitCell3D(convMat, supercell);

supercell(3 : 5, :) = convMat * supercell(3 : 5, :);

WriteLammpsData('tests/lammps/AlGaAs2_3x1x2.lammps', convMat, supercell([1, 3 : 5], :));