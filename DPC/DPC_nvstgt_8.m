% DPC_nvstgt_8.m -- Phase retrieval using MOI
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
% RCoord = sqrt((X-0).^2 + (Y-0).^2);
% % TestPot = CircTable(RCoord, 1, 4, 1500);
% TestPot = CircParabola(RCoord, 1000, 4) + CircParabola(RCoord, 2000, 2);
TestPot = ProjectedPotential(Lx, Ly, Nx, Ny, 14, 0, 0);
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

[OriPhaseGradx, OriPhaseGrady] = gradient(TestPot / 1e3, dx, dy);
EleFilStr = sqrt(OriPhaseGradx.^2 + OriPhaseGrady.^2);
% show PhaseGrad:
figure;
subplot(1, 2, 1);
imagesc(x, y, EleFilStr);
colorbar; axis square;
title('Electric field strength');
subplot(1, 2, 2);
plot(x, EleFilStr(Ny / 2 + 1, : ));
xlabel('x (Angs.)'); ylabel('KeV');
%% Scanning parameter setting:
Scan_Lx = 3;
Scan_Ly = 3;
Scan_Nx = 16;
Scan_Ny = 16;
Scan_dx = Scan_Lx / Scan_Nx;
Scan_dy = Scan_Ly / Scan_Ny;
Scan_x = -Scan_Lx / 2 : Scan_dx : Scan_Lx / 2 -Scan_dx;
Scan_y = -Scan_Ly / 2 : Scan_dy : Scan_Ly / 2 -Scan_dy;
[Scan_X, Scan_Y] = meshgrid(Scan_x, Scan_y);
%% Scanning and data collection:
WB = waitbar(0, 'Scanning initializing...');
TotalCompNum = Scan_Ny * Scan_Nx;
MOI = zeros(Scan_Ny, Scan_Nx);

for ScanY_Idx = 1 : Scan_Ny
    yp = Scan_y(ScanY_Idx);
    for ScanX_Idx = 1 : Scan_Nx
        xp = Scan_x(ScanX_Idx);
        Probe = ProbeCreate(Params, xp, yp, Lx, Ly, Nx, Ny);
        TransWave = Probe .* TransFunc;
        RecoWave = ifftshift(fft2(fftshift(TransWave))) * dx * dy;
        RecoInten = abs(RecoWave.^2);
        MOI(ScanY_Idx, ScanX_Idx) = sum(sum(RecoInten .* Fsqu)) / (Lx * Ly);
        
        CurCompNum = (ScanY_Idx - 1) * Scan_Nx + ScanX_Idx;
        waitbar(roundn(CurCompNum / TotalCompNum, -3), WB, [num2str(roundn(CurCompNum / TotalCompNum, -3) * 100), ' %']);
    end
end
delete(WB);

% Shift MOI:
MOI = MOI - min(MOI( : ));
% Show MOI:
figure;
subplot(1, 2, 1);
imagesc(Scan_x, Scan_y, MOI);
colorbar; axis square;
title('MOI');
subplot(1, 2, 2);
plot(Scan_x, MOI(Scan_Ny / 2 + 1, : ));
xlabel('x (Angs.)');

% save(strcat(SaveDir, 'MOIx.txt'), 'MOIx', '-ascii', '-double', '-tabs');
%% Retrieval using MOI:
ProbeRef = ProbeCreate(Params, 0, 0, Scan_Lx, Scan_Ly, Scan_Nx, Scan_Ny);
ProbeRefInten = abs(ProbeRef.^2);
FT_ProbeRefInten = fft2(fftshift(ProbeRefInten));
FT_ProbeRefInten(FT_ProbeRefInten == 0) = 1e-16;

FT_MOI = fft2(fftshift(MOI));

TempVar_1 = FT_MOI ./ FT_ProbeRefInten;
TempVar_2 = ifftshift(ifft2(TempVar_1)) / (Scan_dx * Scan_dy);
PhaseGrad = 2 * pi * real(sqrt(TempVar_2));

% show PhaseGrad:
figure;
subplot(1, 2, 1);
imagesc(Scan_x, Scan_y, PhaseGrad);
colorbar; axis square;
title('PhaseGrad');
subplot(1, 2, 2);
plot(Scan_x, PhaseGrad(Scan_Ny / 2 + 1, : ));
xlabel('x (Angs.)');

%Compute ref EF:
RefRCoord = sqrt((Scan_X-0).^2 + (Scan_Y-0).^2);
RefPot = CircParabola(RefRCoord, 1000, 4) + CircParabola(RefRCoord, 2000, 2); % in V - Angs;
[RefPotGradx, RefPotGrady] = gradient(RefPot / 1e3, Scan_dx, Scan_dy); % in KeV 
RefEF = sqrt(RefPotGradx.^2 + RefPotGrady.^2);

RetriEF = PhaseGrad / InterCoeff;
% Show RetriEF:
figure;
subplot(1, 2, 1);
imagesc(Scan_x, Scan_y, RetriEF);
colorbar; axis square;
title('RetriEF');
subplot(1, 2, 2);
plot(Scan_x, RetriEF(Scan_Ny / 2 + 1, : ), '-b', Scan_x, RefEF(Scan_Ny / 2 + 1, : ), '--r');
xlabel('x (Angs.)'); ylabel('KeV');
legend('Retri', 'Ori');