% test_MonoProjPot_imtl_X_1.m
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

atomType = 12;
atomNum = 10;
eleProp = ones(1, atomNum);
xyCoords = zeros(2, atomNum);
xyCoords(1, :) = Lx * rand(1, atomNum) - Lx / 2;

projPot = MonoProjPot_imtl_X(atomType, eleProp, xyCoords, Lx, Ly, Nx, Ny);
figure;
imagesc(x, y, projPot);
colormap('gray');
axis square;
xlabel('$x (\AA)$', 'interpreter', 'latex');
ylabel('$y (\AA)$', 'interpreter', 'latex');

figure;
plot(x, projPot(Ny / 2 + 1, :) / 1e3);
hold on;
scatter(xyCoords(1, :), xyCoords(2, :), 'filled');
hold off;
xlabel('$x (\AA)$', 'interpreter', 'latex');
ylabel('$V_z (kV-\AA)$', 'interpreter', 'latex');