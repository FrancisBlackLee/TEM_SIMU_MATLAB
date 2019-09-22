% DPC_nvstgt_7.m -- verify the relation between phase and MOI
clc;
close all;
clear all;
% SaveDir = 'D:\Francis. B. Lee\Practice\Conventional Multislice in MATLAB\DPC\nvstgt_data\MOI_20190914\';
%% Sampling setting:
Lx = 15;
Ly = Lx;
Nx = 1024;
Ny = 1024;
dx = Lx / Nx;
dy = Ly / Ny;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
[X, Y] = meshgrid(x, y);
fx = -1 / (2*dx) : 1 / Lx : 1 / (2*dx) - 1 / Lx;
fy = -1 / (2*dy) : 1 / Ly : 1 / (2*dy) - 1 / Ly;
[FX, FY] = meshgrid(fx, fy);
FsquX = FX.^2;
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
% TestPot = CircTable(RCoord, 1, 4, 1500);
TestPot = CircParabola(RCoord, 2000, 4);
% TestPot = ProjectedPotential(Lx, Ly, Nx, Ny, 14, 0, 0);
TransFunc = exp(1i * InterCoeff * TestPot / 1e3);

% Show TestPot:
figure;
subplot(1, 2, 1);
imagesc(x, y, TestPot);
colormap('gray'); axis square;
title('Projected potential');
subplot(1, 2, 2);
plot(x, TestPot(Ny / 2 + 1, : ) / 1e3);
xlabel('x (Angs.)'); ylabel('Proj.Pot. (KV-Angs.)');

% Phase gradient:
[PhaseGradx, PhaseGrady] = gradient(InterCoeff * TestPot / 1e3, dx, dy);
% show PhaseGradx:
figure;
subplot(1, 2, 1);
imagesc(x, y, PhaseGradx);
colorbar; axis square;
title('x gradient of phase of the transmission function');
subplot(1, 2, 2);
plot(x, PhaseGradx(Ny / 2 + 1, : ));
xlabel('x (Angs.)');

% Prepare Probe:
Probe = ProbeCreate(Params, 0, 0, Lx, Ly, Nx, Ny);
ProbeCG = conj(Probe);
[ProbeGradx, ProbeGrady] = gradient(Probe, dx, dy);
[ProbeCGGradx, ProbeCGGrady] = gradient(ProbeCG, dx, dy);

% save(strcat(SaveDir, 'ProjPot.txt'), 'TestPot', '-ascii', '-double', '-tabs');
% save(strcat(SaveDir, 'PhaseGradx.txt'), 'PhaseGradx', '-ascii', '-double', '-tabs');
%% Solution 1:
subComp1_1 = 3 * 1i * ProbeGradx .* ProbeCG + 1i * Probe .* ProbeCGGradx;
FT_subComp1_1 = fft2(fftshift(subComp1_1)) * dx * dy;
subComp1_2 = PhaseGradx;
FT_subComp1_2 = fft2(fftshift(subComp1_2)) * dx * dy;
subComp1_3 = abs(Probe.^2);
FT_subComp1_3 = fft2(fftshift(subComp1_3)) * dx * dy;
subComp1_4 = PhaseGradx.^2;
FT_subComp1_4 = fft2(fftshift(subComp1_4)) * dx * dy;

Comp1_1 = ifftshift(ifft2(FT_subComp1_1 .* FT_subComp1_2)) / (dx * dy);
Comp1_2 = ifftshift(ifft2(FT_subComp1_3 .* FT_subComp1_4)) / (dx * dy);

MOI_sol_1 = -real(Comp1_1 - Comp1_2) / (4 * pi^2);
% show solution 1:
figure;
subplot(1, 2, 1);
imagesc(x, y, MOI_sol_1);
colorbar; axis square;
title('MOI solution 1');
subplot(1, 2, 2);
plot(x, MOI_sol_1(Ny / 2 + 1, : ));
xlabel('x (Angs.)');

% save(strcat(SaveDir, 'MOI_sol_1.txt'), 'MOI_sol_1', '-ascii', '-double', '-tabs');
%% Solution 2:
subComp2_1 = abs(Probe.^2);
FT_subComp2_1 = fft2(fftshift(subComp2_1)) * dx * dy;
subComp2_2 = PhaseGradx.^2;
FT_subComp2_2 = fft2(fftshift(subComp2_2)) * dx * dy;
subComp2_3 = Probe .* ProbeCGGradx - ProbeGradx .* ProbeCG;
FT_subComp2_3 = fft2(fftshift(subComp2_3)) * dx * dy;
subComp2_4 = PhaseGradx;
FT_subComp2_4 = fft2(fftshift(subComp2_4)) * dx * dy;

Comp2_1 = ifftshift(ifft2(FT_subComp2_1 .* FT_subComp2_2)) / (dx * dy);
Comp2_2 = ifftshift(ifft2(FT_subComp2_3 .* FT_subComp2_4)) / (dx * dy);

MOI_sol_2 = -real(Comp2_1 + 1i * Comp2_2) / (4 * pi^2);
% show solution 2:
figure;
subplot(1, 2, 1);
imagesc(x, y, MOI_sol_2);
colorbar; axis square;
title('MOI solution 2');
subplot(1, 2, 2);
plot(x, MOI_sol_2(Ny / 2 + 1, : ));
xlabel('x (Angs.)');

% save(strcat(SaveDir, 'MOI_sol_2.txt'), 'MOI_sol_2', '-ascii', '-double', '-tabs');
%% Scanning parameter setting:
Scan_Lx = Lx / 1.5;
Scan_Ly = Ly / 1.5;
Scan_Nx = 32;
Scan_Ny = 32;
Scan_dx = Scan_Lx / Scan_Nx;
Scan_dy = Scan_Ly / Scan_Ny;
Scan_X = -Scan_Lx / 2 : Scan_dx : Scan_Lx / 2 -Scan_dx;
Scan_Y = -Scan_Ly / 2 : Scan_dy : Scan_Ly / 2 -Scan_dy;
%% Scanning and data collection:
WB = waitbar(0, 'Scanning initializing...');
TotalCompNum = Scan_Ny * Scan_Nx;
MOIx = zeros(Scan_Ny, Scan_Nx);
for ScanY_Idx = 1 : Scan_Ny
    yp = Scan_Y(ScanY_Idx);
    for ScanX_Idx = 1 : Scan_Nx
        xp = Scan_X(ScanX_Idx);
        Probe = ProbeCreate(Params, xp, yp, Lx, Ly, Nx, Ny);
        TransWave = Probe .* TransFunc;
        RecoWave = ifftshift(fft2(fftshift(TransWave))) * dx * dy;
        RecoInten = abs(RecoWave.^2);
        MOIx(ScanY_Idx, ScanX_Idx) = sum(sum(RecoInten .* FsquX)) / (Lx * Ly);
        CurCompNum = (ScanY_Idx - 1) * Scan_Nx + ScanX_Idx;
        waitbar(roundn(CurCompNum / TotalCompNum, -3), WB, [num2str(roundn(CurCompNum / TotalCompNum, -3) * 100), ' %']);
    end
end
delete(WB);

% Show MOI:
figure;
subplot(1, 2, 1);
imagesc(Scan_X, Scan_Y, MOIx);
colorbar; axis square;
title('MOIx');
subplot(1, 2, 2);
plot(Scan_X, MOIx(Scan_Ny / 2 + 1, : ));
xlabel('x (Angs.)');

% save(strcat(SaveDir, 'MOIx.txt'), 'MOIx', '-ascii', '-double', '-tabs');