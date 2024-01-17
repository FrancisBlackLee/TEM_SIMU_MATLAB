% prototype_triclinic_stem_imaging.m
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

    % figure;
    % subplot(1, 2, 1);
    % mesh(xMesh, yMesh, tmpProjPot);
    % view(0, 90);
    % subplot(1, 2, 2);
    % PlotUnitCell3D(convMat, slices{iSlice});
end

nStack = 10;

%% scanning probe:
[mappedFxMesh, mappedFyMesh] = MsbtFreqMeshTiltedToPlane(fxMesh, fyMesh, wavLen, ...
    invRotMat);
aperture = MeshedCircApert(mappedFxMesh, mappedFyMesh, wavLen, c2Angle);
detector = MeshedAnnularDetector(mappedFxMesh, mappedFyMesh, wavLen, 76, 200);

% scanning points
scanN1 = round(norm(a1(1 : 2)) / 0.2) + 1;
scanN2 = round(norm(a2(1 : 2)) / 0.2) + 1;
scanFracA1 = linspace(0, 1, scanN1) - 0.5;
scanFracA2 = linspace(0, 1, scanN2) - 0.5;
[scanFracA1Mesh, scanFracA2Mesh] = meshgrid(scanFracA1, scanFracA2);
scanX = scanFracA1Mesh * a1(1) + scanFracA2Mesh * a2(1);
scanY = scanFracA1Mesh * a1(2) + scanFracA2Mesh * a2(2);
% figure('units','normalized','outerposition',[0 0 1 1]);
haadf = zeros(scanN2, scanN1);
wb = waitbar(0, 'init scanning...');
for scanI2 = 1 : scanN2
    parfor scanI1 = 1 : scanN1
        px = scanX(scanI2, scanI1);
        py = scanY(scanI2, scanI1);
        probe = MeshedProbe(aberrs, wavLen, aperture, px, py, mappedFxMesh, ...
            mappedFyMesh);
        % probeI = abs(probe.^2);
        % 
        % mesh(xMesh, yMesh, probeI);
        % xlabel('x (\AA)', 'Interpreter', 'latex');
        % ylabel('y (\AA)', 'Interpreter', 'latex');
        % title(['probe on tilted plane(x = ', num2str(px), ', y = ', num2str(py), ')']);
        % view(0, 90);
        % drawnow;
        % pause(0.1);

        wave = MsbtMultislice(probe, wavLen, fxMesh, fyMesh, a3p, transFuncs, sliceDists, nStack);
        wave = ifftshift(fft2(fftshift(wave)));
        waveI = abs(wave.^2);
        haadf(scanI2, scanI1) = VtemlabPacbedQstem(waveI, detector);
    end

    waitbar(scanI2 / scanN2, wb, [num2str(scanI2), '/', num2str(scanN2)]);
end
close(wb);

figure;
mesh(scanX, scanY, haadf, 'FaceColor', 'flat');
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('HAADF');
axis equal;
view(0, 90);

%% repeat the haadf of unit cell:
repA1 = 10;
repA2 = 10;
repHaadf = repmat(haadf(1 : scanN2 - 1, 1 : scanN1 - 1), repA2, repA1);
repScanFracA1 = InitAxis(repA1, repA1 * (scanN1 - 1));
repScanFracA2 = InitAxis(repA2, repA2 * (scanN2 - 1));
[repScanFracA1Mesh, repScanFracA2Mesh] = meshgrid(repScanFracA1, repScanFracA2);
repScanX = repScanFracA1Mesh * a1(1) + repScanFracA2Mesh * a2(1);
repScanY = repScanFracA1Mesh * a1(2) + repScanFracA2Mesh * a2(2);

figure;
mesh(repScanX, repScanY, repHaadf, 'FaceColor', 'flat');
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
title('Repeated HAADF');
axis equal;
view(0, 90);