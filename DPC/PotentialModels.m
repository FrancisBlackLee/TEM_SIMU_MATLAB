clc;
clear all;
close all;
%% Compare the models of atomic potential
% Sampling parameters:
Lx = 2;
Ly = Lx;
Nx = 1024;
Ny = Nx;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
% Model in Advanced Computing in Electron Microscopy
Vstd = ProjectedPotential(Lx, Ly, Nx, Ny, 10, 0, 0);
% Exponential model:
Vmax = max(max(Vstd));
R = sqrt(X.^2 + Y.^2);
Vexp = Vmax * exp(-log(Vmax) * R.^0.75);
% Plot the profiles of these models:
figure;
plot(x, Vstd(Ny / 2 + 1, : ), x, Vexp(Ny / 2 + 1, : ));
legend('std', 'exp');