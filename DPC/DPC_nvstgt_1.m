% DPC_nvstgt_1: the form of linear field of potential
% Silicon [1 1 0] is used
clc;
close all;
clear all;
%% Prepare the specimen
Lattice_Const = [3.84, 5.43]; % [a b]
SliceDist = [1.9198, 1.9198]; % distance between each slice
Cell_Num = [6, 4]; % expand the unit cell by Expan_Nx = 3 and Expan_Ny = 2, adaptive
DistError = 1e-2;
% Laters: Each column for an atom
SliceA = [0, 0.5; 0, 0.75];
SliceB = [0, 0.5; 0.25, 0.5];
% Expansion
SliceA = SquareLattExpan(SliceA, Lattice_Const, Cell_Num);
SliceB = SquareLattExpan(SliceB, Lattice_Const, Cell_Num);
%% basic settings
% sampling:
Lx = Cell_Num(1) * Lattice_Const(1);
Ly = Cell_Num(2) * Lattice_Const(2);
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
R = sqrt(X.^2 + Y.^2);
%% Slice potentials test
% Layer A:
ProjPotA = ProjectedPotential(Lx, Ly, Nx, Ny, 14 * ones(size(SliceA, 2), 1), SliceA(1, :), SliceA(2, :));
% Layer B:
ProjPotB = ProjectedPotential(Lx, Ly, Nx, Ny, 14 * ones(size(SliceB, 2), 1), SliceB(1, :), SliceB(2, :));
% Gradients of the potential fields
[GradPotAx, GradPotAy] = gradient(ProjPotA, dx, dy);
[GradPotBx, GradPotBy] = gradient(ProjPotB, dx, dy);
% scalar product of potential gradient field and position vector field
TestPotA = R .* (GradPotAx + GradPotAy);
TestPotB = R .* (GradPotBx + GradPotBy);
% Comparison
% PotA:
figure;
subplot(2, 2, 1);
imagesc(x, y, ProjPotA);
colormap('gray');
axis square;
title('ProjPotA');
subplot(2, 2, 2);
imagesc(x, y, TestPotA);
colormap('gray');
axis square;
title('TestPotA');
subplot(2, 2, 3);
plot(x, ProjPotA(Ny / 2 + 1, : ));
title('ProjPotA');
subplot(2, 2, 4);
plot(x, TestPotA(Ny / 2 + 1, : ));
title('TestPotA');
% PotB:
figure;
subplot(2, 2, 1);
imagesc(x, y, ProjPotB);
colormap('gray');
axis square;
title('ProjPotB');
subplot(2, 2, 2);
imagesc(x, y, TestPotB);
colormap('gray');
axis square;
title('TestPotB');
subplot(2, 2, 3);
plot(y, ProjPotB( : , Nx / 2 + 1));
title('ProjPotB');
subplot(2, 2, 4);
plot(y, TestPotB( : , Nx / 2 + 1));
title('TestPotB');