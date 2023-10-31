% test_MeshedProjPotKs.m
clc;
clear;
close all;
%% hexagonal unit cell:
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

%% orthogonal unit cell:
orthoUnitCell = [6, 1, 0, 0, 1/2;...
    6, 1, 1/2, 1/6, 1/2;...
    6, 1, 1/2, 1/2, 1/2;...
    6, 1, 0, 2/3, 1/2]';
orthoLattConsts = a * [1, sqrt(3), 1];
orthoConvMat = zeros(3);
orthoConvMat(1, 1) = orthoLattConsts(1);
orthoConvMat(2, 2) = orthoLattConsts(2);
orthoConvMat(3, 3) = orthoLattConsts(3);

expanNum = [3, 2];
lx = expanNum(1) * orthoLattConsts(1);
ly = expanNum(2) * orthoLattConsts(2);
nx = 512;
ny = 512;
x = InitAxis(lx, nx);
y = InitAxis(ly, ny);

figure;
PlotUnitCell2D(orthoConvMat, orthoUnitCell);

orthoProjPot = MultiProjPot_conv_X(orthoUnitCell, [3, 2], orthoLattConsts, lx, ly, nx, ny);

figure;
mesh(x, y, orthoProjPot);
xlabel('x (\AA)', 'Interpreter', 'latex');
ylabel('y (\AA)', 'Interpreter', 'latex');
zlabel('$Projected potential (V \cdot \AA)$', 'Interpreter', 'latex');