% prototype_triclinic_plane_wave.m
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
sliceDists = norm(a3p) * sliceDists;

%% meshes:
n1 = 1024;
n2 = 1024;
[xMesh, yMesh, fxMesh, fyMesh] = InitNonorthoMesh2D(a1p, a2p, na1, na2, n1, n2);

%% tem settings:
keV = 300;
wavLen = HighEnergyWavLen_X(keV);

pol = acos(dot(a3p, [0,0,1]') / norm(a3p));
azi = atan2(a3p(2), a3p(1));

fx0 = pol * cos(azi) / wavLen;
fy0 = pol * sin(azi) / wavLen;

inWave = exp(1i * 2 * pi * (fx0 * xMesh + fy0 * yMesh));

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

%% propagation kernels:
ku = a3p / norm(a3p);
propKernels = 1i * ones(n2, n1, nSlice);
for iSlice = 1 : nSlice
    propKernels(:, :, iSlice) = MsbtPropKernel(fxMesh, fyMesh, ku, wavLen, sliceDists(iSlice));
end

%% multilsice:
nStack = 20;
wave = inWave;
figure;
for iStack = 1 : nStack
    for iSlice = 1 : nSlice
        wave = wave .* transFuncs(:, :, iSlice);
        wave = ifftshift(ifft2(fftshift(propKernels(:, :, iSlice)) .*...
            fft2(fftshift(wave))));
        
        waveI = abs(wave.^2);
        mesh(xMesh, yMesh, waveI);
        view(0, 90);
        drawnow;
    end
end