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
colormap('jet');
axis equal;
axis tight;
title('Nonorthogonal mesh: CBED');
view(0, 90);

%% load ref CBED:
refLy = 38.7038 * 2;
refLx = 38.3069 * 2;
refNx = 512 * 2;
refNy = 512 * 2;
refFx = InitFreqAxis(refLx, refNx);
refFy = InitFreqAxis(refLy, refNy);
refCbed = ReadBinaryFile('tests/nonorthogonal/MoS2_110_cbed.bin', [refNy, refNx], 'row', 'float');
refCbed = refCbed / sum(refCbed, 'all');

figure;
imagesc(refFy, refFx, refCbed');
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
colormap('jet');
axis equal;
axis tight;
title('Orthogonal mesh: CBED');
view(0, 90);