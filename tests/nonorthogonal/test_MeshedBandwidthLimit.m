% test_MeshedBandwidthLimit.m
clc;
clear;
close all;
%% projected potential:
a = 2;
a1 = a * [1/2, sqrt(3)/2, 0.0];
a2 = a * [1/2, -sqrt(3)/2, 0.0];
a3 = [0, 0, a];
hexUnitCell = [6, 1, 1/3, 2/3, 1/2;...
    6, 1, 2/3, 1/3, 1/2]';
hexConvMat = [a1', a2', a3'];

na1 = 3;
na2 = 3;
tiles = [na1, na2, 1];
hexSuperCell = TileUnitCell(hexUnitCell, tiles);

figure;
PlotUnitCell2D(hexConvMat, hexSuperCell);

n1 = 512;
n2 = 512;
[xMesh, yMesh, fxMesh, fyMesh] = InitNonorthoMesh2D(a1, a2, na1, na2, n1, n2);
stretchCoeff = ScattFacStretchCoeff(a1, a2, na1, na2, n1, n2);

hexSuperCell(3, :) = hexSuperCell(3, :) - na1 / 2;
hexSuperCell(4, :) = hexSuperCell(4, :) - na2 / 2;
hexSuperCell(3 : 5, :) = hexConvMat * hexSuperCell(3 : 5, :);

hexProjPot = MeshedProjPotKs(fxMesh, fyMesh, stretchCoeff, hexSuperCell(1, :), ...
    hexSuperCell(2, :), hexSuperCell(3, :), hexSuperCell(4, :));

figure;
mesh(xMesh, yMesh, hexProjPot);
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
zlabel('$Projected potential (V \cdot \AA)$', 'Interpreter', 'latex');

%% bandwidth limit transmission function:
keV = 300;
interCoeff = InteractionCoefficient(keV);
hexTransFunc = exp(1i * interCoeff * hexProjPot / 1.0e3);

figure;
mesh(xMesh, yMesh, angle(hexTransFunc));
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
zlabel('$TF phase (V \cdot \AA)$', 'Interpreter', 'latex');

bwl = NonorthoMeshBwlFreq(a1, a2, na1, na2, n1, n2, 2/3);
bwlHexTransFunc = MeshedBandwidthLimit(hexTransFunc, fxMesh, fyMesh, bwl);

figure;
mesh(xMesh, yMesh, angle(bwlHexTransFunc));
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
zlabel('$Bandwidth limited TF phase (V \cdot \AA)$', 'Interpreter', 'latex');