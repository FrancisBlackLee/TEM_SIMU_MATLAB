% prototype_triclinic_stem_multislice.m
clc;
clear;
close all;
%% load unit cell:
crysInfo = LoadCif('tests/nonorthogonal/Gd2B4O9.cif');
[cellLengths, cellAngles] = ExtractCellConstFromCrysInfo(crysInfo);
[unitCell, convMat, cutoff] = CrystalAdvisor('tests/nonorthogonal/Gd2B4O9.cif', ...
    [0, 0, 1], true);

figure;
PlotUnitCell3D(convMat, unitCell);
title('unit cell');

[a1, a2, a3] = ConvMatToBases(convMat);
[a1p, a2p, rotMat, invRotMat] = NonorthoBasesReorient(a1, a2);
convMatP = [reshape(a1p, [], 1), reshape(a2p, [], 1), zeros(3, 1)];

a3p = rotMat * reshape(a3, [], 1);
a3Proj = dot(a1, cross(a2, a3)) / norm(cross(a1, a2));

%% build a supercell:
na1 = 3;
na2 = 3;
supercell = TileUnitCell(unitCell, [na1, na2, 1]);

figure;
PlotUnitCell3D(convMat, supercell);
title('supercell (tiled on a1-a2 plane)');

%% slice the supercell:
[slices, sliceDists, extraSlice] = CrystalSlicing_X(supercell, supercell, ...
    0.3, 1, 0, 0);
nSlice = length(sliceDists);
sliceDists = a3Proj * sliceDists;

%% STEM settings:
keV = 300;
wavLen = HighEnergyWavLen_X(keV);
aberrs = InitObjectiveLensAberrations_X();
c2Angle = 20;

%% meshes:
n1 = 1024;
n2 = 1024;
[xMesh, yMesh, fxMesh, fyMesh] = InitNonorthoMesh2D(a1p, a2p, na1, na2, n1, n2);

%% transmission function:
interCoeff = InteractionCoefficient(keV);
stretchCoeff = ScattFacStretchCoeff(a1p, a2p, na1, na2, n1, n2);
transFuncs = 1i * ones(n2, n1, nSlice);
bwl = NonorthoMeshBwlFreq(a1p, a2p, na1, na2, n1, n2, 2/3);

for iSlice = 1 : nSlice
    tmpSlice = slices{iSlice};
    tmpSlice(3 : 4, :) = tmpSlice(3 : 4, :) - [na1; na2] / 2;
    tmpSlice(3 : 5, :) = convMatP * tmpSlice(3 : 5, :);
    tmpProjPot = MeshedProjPotKs(fxMesh, fyMesh, stretchCoeff, tmpSlice(1, :), ...
        tmpSlice(2, :), tmpSlice(3, :), tmpSlice(4, :));
    transFuncs(:, :, iSlice) = exp(1i * interCoeff * tmpProjPot / 1.0e3);
    transFuncs(:, :, iSlice) = MeshedBandwidthLimit(transFuncs(:, :, iSlice), ...
        fxMesh, fyMesh, bwl);

%     figure;
%     subplot(1, 2, 1);
%     mesh(xMesh, yMesh, tmpProjPot);
%     subplot(1, 2, 2);
%     PlotUnitCell3D(convMat, slices{iSlice});
end

%% probe:
[mappedFxMesh, mappedFyMesh] = MsbtFreqMeshTiltedToPlane(fxMesh, fyMesh, wavLen, ...
    invRotMat);
aperture = MeshedCircApert(mappedFxMesh, mappedFyMesh, wavLen, c2Angle);
probe = MeshedProbe(aberrs, wavLen, aperture, 0, 0, mappedFxMesh, mappedFyMesh);
probeI = abs(probe.^2);

figure;
mesh(xMesh, yMesh, probeI);
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('probe on tilted plane');
view(0, 90);

%% propagation kernels:
ku = a3p / norm(a3p);
propKernels = 1i * ones(n2, n1, nSlice);
for iSlice = 1 : nSlice
    propKernels(:, :, iSlice) = MsbtPropKernel(fxMesh, fyMesh, ku, wavLen, sliceDists(iSlice));
end

%% multilsice:
nStack = 20;
wave = probe;
for iStack = 1 : nStack
    for iSlice = 1 : nSlice
        wave = wave .* transFuncs(:, :, iSlice);
        wave = ifftshift(ifft2(fftshift(propKernels(:, :, iSlice)) .*...
            fft2(fftshift(wave))));
        
%         waveI = abs(wave.^2);
%         mesh(xMesh, yMesh, waveI);
%         view(0, 90);
    end
end

%% CBED on tilted plane:
wave = ifftshift(fft2(fftshift(wave)));
waveI = abs(wave.^2);
logWaveI = log(1 + 1000 * waveI / max(waveI, [], 'all'));
figure;
mesh(mappedFxMesh, mappedFyMesh, logWaveI);
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
title('CBED on tilted plane (log scale)');
axis equal;
view(0, 90);

%% transform wave from tiled plane to detector plane:
waveReal = real(wave);
waveImag = imag(wave);

interpReal = scatteredInterpolant(fxMesh(:), fyMesh(:), waveReal(:), 'linear', 'nearest');
interpImag = scatteredInterpolant(fxMesh(:), fyMesh(:), waveImag(:), 'linear', 'nearest');

[detFxMesh, detFyMesh] = MsbtFreqMeshTiltedToPlane(fxMesh, fyMesh, wavLen, rotMat);
detWaveReal = interpReal(detFxMesh(:), detFyMesh(:));
detWaveImag = interpImag(detFxMesh(:), detFyMesh(:));

detWaveReal = reshape(detWaveReal, n2, n1);
detWaveImag = reshape(detWaveImag, n2, n1);

detWave = detWaveReal + 1i * detWaveImag;
detWaveI = abs(detWave.^2);

figure;
mesh(fxMesh, fyMesh, detWaveI);
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
title('CBED on detector plane');
axis equal;
view(0, 90);

logDetWaveI = log(1 + 1000 * detWaveI / max(detWaveI, [], 'all'));
figure;
mesh(fxMesh, fyMesh, logDetWaveI);
xlabel('$f_x$ (1 / \AA)', 'Interpreter', 'latex');
ylabel('$f_y$ (1 / \AA)', 'Interpreter', 'latex');
title('CBED on detector plane (log scale)');
axis equal;
view(0, 90);

%% why not map detector from detector plane to tilted plane:
detector = MeshedAnnularDetector(mappedFxMesh, mappedFyMesh, wavLen, 76, 200);

figure;
mesh(mappedFxMesh * wavLen * 1e3, mappedFyMesh * wavLen * 1e3, detector);
xlabel('$\beta_x$ (mrad)', 'Interpreter', 'latex');
ylabel('$\beta_y$ (mrad)', 'Interpreter', 'latex');
title('HAADF detector on tilted plane (log scale)');
axis equal;
view(0, 90);