% test_MeshedMultislice.m
clc;
clear;
close all;
%% hex unit cell and supercell:
crysInfo = LoadCif('tests/nonorthogonal/MoS2.cif');
[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
hexUnitCell = ExtractAtomSiteFromCrysInfo(crysInfo);
hexConvMat = ConversionMatrix(cellLengths, cellAngles);

figure;
PlotUnitCell3D(hexConvMat, hexUnitCell);

a1 = hexConvMat(:, 1)';
a2 = hexConvMat(:, 2)';
a3 = hexConvMat(:, 3)';

na1 = 16;
na2 = 16;

n1 = 512;
n2 = 512;

hexSuperCell = TileUnitCell(hexUnitCell, [na1, na2, 1]);
hexSuperCell(3, :) = hexSuperCell(3, :) - na1 / 2;
hexSuperCell(4, :) = hexSuperCell(4, :) - na2 / 2;
hexSuperCell(3 : 5, :) = hexConvMat * hexSuperCell(3 : 5, :);
[hexSlices, hexSliceDists, extraSlice] = CrystalSlicing_X(hexSuperCell, hexSuperCell,...
    1.5, norm(a3), 0, 0);
nHexSlice = length(hexSliceDists);

%% transmission functions:
keV = 300;
interCoeff = InteractionCoefficient(keV);
[xMesh, yMesh, fxMesh, fyMesh] = InitNonorthoMesh2D(a1, a2, na1, na2, n1, n2);
stretchCoeff = ScattFacStretchCoeff(a1, a2, na1, na2, n1, n2);
hexTransFuncs = 1i * ones(n1, n2, nHexSlice);
bwl = NonorthoMeshBwlFreq(a1, a2, na1, na2, n1, n2, 2/3);

for iHexSlice = 1 : nHexSlice
    tmpProjPot = MeshedProjPotKs(fxMesh, fyMesh, stretchCoeff, ...
        hexSlices{iHexSlice}(1, :), hexSlices{iHexSlice}(2, :), ...
        hexSlices{iHexSlice}(3, :), hexSlices{iHexSlice}(4, :));
    hexTransFuncs(:, :, iHexSlice) = exp(1i * interCoeff * tmpProjPot / 1.0e3);
    hexTransFuncs(:, :, iHexSlice) = MeshedBandwidthLimit(hexTransFuncs(:, :, iHexSlice), ...
        fxMesh, fyMesh, bwl);
end

%% incident wave:
numApert = 3.0;
wavLen = HighEnergyWavLen_X(keV);
aberrs = InitObjectiveLensAberrations_X();
aperture = MeshedCircApert(fxMesh, fyMesh, wavLen, numApert);
probe = MeshedProbe(aberrs, wavLen, aperture, 0.0, 0.0, fxMesh, fyMesh);

%% CBED:
nStack = 37;
wave = MeshedMultislice(probe, wavLen, fxMesh, fyMesh, hexTransFuncs, hexSliceDists, nStack);
wave = ifftshift(fft2(fftshift(wave)));
waveI = abs(wave.^2);
waveI = waveI / sum(waveI, 'all');
logWaveI = log(1 + 100 * waveI / max(waveI, [], 'all'));

figure;
mesh(fxMesh, fyMesh, waveI);
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
axis equal;
axis tight;
title('Nonorthogonal mesh: CBED');

%% load ref CBED:
lx = 38.7038;
ly = 38.3069;
nx = 512;
ny = 512;
x = InitAxis(lx, nx);
y = InitAxis(ly, ny);
refCbed = ReadBinaryFile('tests/nonorthogonal/MoS2_110_cbed.bin', [ny, nx], 'row', 'float');
refCbed = refCbed / sum(refCbed, 'all');

figure;
imagesc(x, y, refCbed);
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
axis equal;
axis tight;
title('Orthogonal mesh: CBED');