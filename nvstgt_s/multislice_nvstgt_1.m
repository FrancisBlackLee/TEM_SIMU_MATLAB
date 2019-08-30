% multislice_nvstgt_1.m test the projected potential generated using
% convolution
% Si [110] sample:
clc;
close all;
clear all;
%% Main
LattConst = [3.84, 5.43]; % [a b]
M = 2;
CellNum = M * [3, 2]; % expand the unit cell by Expan_Nx = 3 and Expan_Ny = 2, adaptive
LayerA = [0, 0.5; 0, 0.75];
LayerB = [0, 0.5; 0.25, 0.5];

Lx = CellNum(1) * LattConst(1);
Ly = CellNum(2) * LattConst(2);
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
ProjPotA = MonoProjPot_conv_0(14, LayerA, CellNum, LattConst, Lx, Ly, Nx, Ny);
ProjPotB = MonoProjPot_conv_0(14, LayerB, CellNum, LattConst, Lx, Ly, Nx, Ny);

% Show the results:
figure;
imagesc(x, y, ProjPotA);
colormap('gray'); axis square;
title('Slice A Proj.Pot.');
figure;
imagesc(x, y, ProjPotB);
colormap('gray'); axis square;
title('Slice B Proj.Pot.');

% Lattice_Const = [3.84, 5.43]; % [a b]
% LayerDist = [1.9198, 1.9198]; % distance between each slice
% M = 10;
% Cell_Num = M * [3, 2]; % expand the unit cell by Expan_Nx = 3 and Expan_Ny = 2, adaptive
% DistError = 1e-2;
% % Laters: Each column for an atom
% LayerA = [0, 0.5; 0, 0.75];
% LayerB = [0, 0.5; 0.25, 0.5];
% % Expansion
% LayerA = SquareLattExpan(LayerA, Lattice_Const, Cell_Num);
% LayerB = SquareLattExpan(LayerB, Lattice_Const, Cell_Num);
% % sampling:
% Lx = Cell_Num(1) * Lattice_Const(1);
% Ly = Cell_Num(2) * Lattice_Const(2);
% Nx = 1024;
% Ny = 1024;
% dx = Lx / Nx;
% dy = Ly / Ny;
% x = -Lx / 2 : dx : Lx / 2 - dx;
% y = -Ly / 2 : dy : Ly / 2 - dy;
% % Layer A:
% Proj_PotA = ProjectedPotential(Lx, Ly, Nx, Ny, 14 * ones(size(LayerA, 2), 1), LayerA(1, :), LayerA(2, :));
% % Layer B:
% Proj_PotB = ProjectedPotential(Lx, Ly, Nx, Ny, 14 * ones(size(LayerB, 2), 1), LayerB(1, :), LayerB(2, :));
% % test
% figure;
% imagesc(x, y, Proj_PotA);
% colormap('gray'); axis square;
% title('slice A Proj.Pot.');
% figure;
% imagesc(x, y, Proj_PotB);
% colormap('gray'); axis square;
% title('slice B Proj.Pot.');