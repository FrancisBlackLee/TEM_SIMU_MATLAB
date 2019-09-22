% DPC_nvstgt_6.m -- the final shot for the old method:
clc;
clear all;
close all;
%% Sampling setting:
Lx = 15;
Ly = Lx;
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
fx = -1 / (2*dx) : 1 / Lx : 1 / (2*dx) - 1 / Lx;
fy = -1 / (2*dy) : 1 / Ly : 1 / (2*dy) - 1 / Ly;
[FX, FY] = meshgrid(fx, fy);
Fsqu = FX.^2 + FY.^2;
%% STEM parameter setting:
Params.KeV = 300;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber
Params.amax = 23;
Params.Cs = 0;
Params.df = 0;
%% Sample preparation:
RCoord = sqrt((X-0).^2 + (Y-0).^2);
% TestPot = CircTable(RCoord, 0, 2, 1500);
TestPot = CircParabola(RCoord, 1000, 2);
% TestPot = ProjectedPotential(Lx, Ly, Nx, Ny, 14, 0, 0);
TransFunc = exp(1i * InterCoeff * TestPot / 1e3);

% Show TestPot:
figure;
subplot(1, 2, 1);
imagesc(x, y, TestPot);
colormap('gray'); axis square;
subplot(1, 2, 2);
plot(x, TestPot(Ny / 2 + 1, : ) / 1e3);
xlabel('x (Angs.)'); ylabel('Proj.Pot. (KV-Angs.)');

% Add Laplacian (Nabla twice) to the projected potential
[PotGradx, PotGrady] = gradient(TestPot / 1e3, dx, dy);
[PotGradxx, PotGradxy] = gradient(PotGradx, dx, dy);
[PotGradyx, PotGradyy] = gradient(PotGrady, dx, dy);
EleField = sqrt(PotGradx.^2 + PotGrady.^2);
PotLap = PotGradxx + PotGradyy;
% Show EleField:
figure;
subplot(1, 2, 1);
imagesc(x, y, EleField);
colormap('gray'); axis square;
title('Electric filed strength');
subplot(1, 2, 2);
plot(x, EleField(Ny / 2 + 1, : ));
xlabel('x (Angs.)'); ylabel('Proj.Ele.F.Str. (KV)');
% Show PotLap:
figure;
subplot(1, 2, 1);
imagesc(x, y, PotLap);
colormap('gray'); axis square;
title('\nabla^{2}V_{z}');
subplot(1, 2, 2);
plot(x, PotLap(Ny / 2 + 1, : ));
xlabel('x (Angs.)'); ylabel('Proj.Ele.Den. (KV / Angs.)');
%% Scanning parameter setting:
Scan_Lx = Lx / 1.5;
Scan_Ly = Ly / 1.5;
Scan_Nx = 128;
Scan_Ny = 128;
Scan_dx = Scan_Lx / Scan_Nx;
Scan_dy = Scan_Ly / Scan_Ny;
Scan_X = -Scan_Lx / 2 : Scan_dx : Scan_Lx / 2 -Scan_dx;
Scan_Y = -Scan_Ly / 2 : Scan_dy : Scan_Ly / 2 -Scan_dy;
%% Scanning and data collection:
WB = waitbar(0, 'Initializing...');
TotalCompNum = Scan_Ny * Scan_Nx;
COMx = zeros(Scan_Ny, Scan_Nx);
COMy = zeros(Scan_Ny, Scan_Nx);
MOI = zeros(Scan_Ny, Scan_Nx);
for ScanY_Idx = 1 : Scan_Ny
    yp = Scan_Y(ScanY_Idx);
    for ScanX_Idx = 1 : Scan_Nx
        xp = Scan_X(ScanX_Idx);
        Probe = ProbeCreate(Params, xp, yp, Lx, Ly, Nx, Ny);
        TransWave = Probe .* TransFunc;
        RecoWave = ifftshift(fft2(fftshift(TransWave))) * dx * dy;
        RecoInten = abs(RecoWave.^2);
        COMx(ScanY_Idx, ScanX_Idx) = sum(sum(RecoInten .* FX)) / (Lx * Ly);
        COMy(ScanY_Idx, ScanX_Idx) = sum(sum(RecoInten .* FY)) / (Lx * Ly);
        MOI(ScanY_Idx, ScanX_Idx) = sum(sum(RecoInten .* Fsqu)) / (Lx * Ly);
        CurCompNum = (ScanY_Idx - 1) * Scan_Nx + ScanX_Idx;
        waitbar(roundn(CurCompNum / TotalCompNum, -3), WB, [num2str(roundn(CurCompNum / TotalCompNum, -3) * 100), ' %']);
    end
end
delete(WB);
%% Retrieval:
[COMxGradx, COMxGrady] = gradient(COMx, 1 / Scan_Lx, 1 / Scan_Ly);
[COMyGradx, COMyGrady] = gradient(COMy, 1 / Scan_Lx, 1 / Scan_Ly);
dCOM = COMxGradx + COMyGrady;

% Show dCOM:
figure;
subplot(1, 2, 1);
imagesc(x, y, dCOM);
colormap('gray'); axis square;
title('dCOM');
subplot(1, 2, 2);
plot(Scan_X, 2 * pi * dCOM(Scan_Ny / 2 + 1, : ) / InterCoeff);
title('2\pi dCOM / \sigma');

%% Retrieval with dCOM:
RetriPot = iLaplacian_0(2 * pi * dCOM, Scan_Lx, Scan_Ly) / InterCoeff;
RetriPot = RetriPot - min(RetriPot(:));
figure;
subplot(1, 2, 1);
imagesc(Scan_X, Scan_Y, RetriPot);
colormap('gray'); axis square;
title('Retrieved potential (KV-Angs.)');
subplot(1, 2, 2);
plot(Scan_X, RetriPot(Scan_Ny / 2 + 1, : ));