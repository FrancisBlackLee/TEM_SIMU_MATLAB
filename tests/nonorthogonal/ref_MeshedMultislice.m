% ref_MeshedMultislice.m
clc;
clear;
close all;
%% unit cell:
[unitCell, convMat, cutoff] = CrystalAdvisor('tests/nonorthogonal/MoS2.cif', [0, 0, 1], true);
unitCell = RemoveSymmetricAtoms(unitCell);

figure;
PlotUnitCell3D(convMat, unitCell);

% write ejkxyz
nAtom = size(unitCell, 2);
unitCell(3 : 5, :) = convMat * unitCell(3 : 5, :);
lattConsts = [norm(convMat(:, 1)), norm(convMat(:, 2)), norm(convMat(:, 3))];
wobbles = zeros(1, nAtom);

WriteEjkXyz('tests/nonorthogonal/MoS2_110.ejkxyz', lattConsts, unitCell, wobbles);

%% cbed:
task = InitCbedGpuTask();

task.specimen.filename = 'tests/nonorthogonal/MoS2_110.ejkxyz';
task.specimen.tile = [7, 12];
task.specimen.nConfig = 1;
task.specimen.depths = 495;

task.beam.keV = 300;

task.nPixel.nx = 512;
task.nPixel.ny = 512;

task.aperture.r = 3.0;

task.probePosition.x = 0.0;
task.probePosition.y = 0.0;

task.output.name = 'tests/nonorthogonal/MoS2_110_cbed.bin';

ExecuteCbedGpuTask(task, 'float');