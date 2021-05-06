% test_MultiProjPot_imtl_X_1.m
clc;
clear;
close all;
%% main:
Lx = 20;
Ly = Lx;
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;

atomNum = 10;
atomTypeCoords = zeros(4, atomNum);
atomTypeCoords(1, :) = 10 + randi(5, 1, atomNum);
atomTypeCoords(2, :) = 1;
atomTypeCoords(3 : 4, :) = Lx * rand(2, atomNum) - Lx / 2;

projPot = MultiProjPot_imtl_X(atomTypeCoords, Lx, Ly, Nx, Ny);
figure;
imagesc(x, y, projPot);
colormap('gray');
hold on;
scatter(atomTypeCoords(3, :), atomTypeCoords(4, :), '+');
axis square;
xlabel('$x (\AA)$', 'interpreter', 'latex');
ylabel('$y (\AA)$', 'interpreter', 'latex');