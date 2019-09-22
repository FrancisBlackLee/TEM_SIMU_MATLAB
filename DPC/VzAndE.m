% Projected Atomic Potential and Electric Field Strength vs Sampling
% interval
clc;
clear all;
close all;
%% Data preparation
Lx = 4;
Ly = Lx;
Nx = 512;
Ny = Nx;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
GoldPot = ProjectedPotential(Lx, Ly, Nx, Ny, 79, 0, 0);
[Ex, Ey] = gradient(GoldPot, dx, dy);
Ex = -Ex; Ey = -Ey;
E = sqrt(Ex.^2 + Ey.^2);
[Exx, Exy] = gradient(Ex, dx, dy);
[Eyx, Eyy] = gradient(Ey, dx, dy);
ChargeDensity = Exx + Eyy;
% Show the results:
subplot(3, 2, 1);
imagesc(x, y, GoldPot);
colormap('gray');
axis square;
title('projected potential of gold atom');
subplot(3, 2, 2);
plot(x, GoldPot(Ny / 2 + 1, : ));
title('projected potential of gold atom');
subplot(3, 2, 3);
imagesc(x, y, E);
colormap('gray');
axis square;
title('electric field strength');
subplot(3, 2, 4);
plot(x, E(Ny / 2 + 1, : ));
title('electric field strength');
subplot(3, 2, 5);
imagesc(x, y, ChargeDensity);
colormap('gray');
axis square;
title('charge density');
subplot(3, 2, 6);
plot(x, ChargeDensity(Ny / 2 + 1, : ));
title('charge density');