% DPC_nvstgt_2: electron probe located at the center of a cone-shape
% potential distribution.
% Image CBED pattern
clc;
close all;
clear all;
%% basic settings
% sampling:
Lx = 10; % in Angstrom
Ly = Lx;
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
% STEM settings:
Params.KeV = 100;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber
Params.amax = 11.43;
Params.Cs = 1.3;
Params.df = 850;
%% Sample preparation
ConeR = sqrt(X.^2 + Y.^2); % in Angstrom
BottomR = 2; % in Angstrom
ConeHeight = 1000; % in Volt-Angstrom
ConePot = -ConeHeight / BottomR * ConeR + ConeHeight;
ConePot(ConePot < 0) =0;
ConeTF = exp(1i * InterCoeff * ConePot / 1e3);
%% Scanning module
Probe = ProbeCreate(Params, 0, 0, Lx, Ly, Nx, Ny);
Trans_Wave = Probe .* ConeTF;
Trans_Wave_Far = ifftshift(fft2(fftshift(Trans_Wave)) * dx * dy);
% DetectInten = log(abs(Trans_Wave_Far.^2));
DetectInten = abs(Trans_Wave_Far.^2);

ImgXStart = (Nx - 128) / 2 + 1;
ImgXEnd = (Nx + 128) / 2;
ImgYStart = (Ny - 128) / 2 + 1;
ImgYEnd = (Ny + 128) / 2;
% Show the detected image:
figure;
imagesc(x(ImgXStart : ImgXEnd), y(ImgYStart : ImgYEnd), DetectInten(ImgYStart : ImgYEnd, ImgXStart : ImgXEnd));
colormap('gray');
axis square;