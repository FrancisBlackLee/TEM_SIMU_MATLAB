% DPC_nvstgt_3: electron probe located at the center of a cone-shape
% potential distribution.
% Image CBED pattern
% Calculate the equivalent rotational inertia
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
Params.KeV = 300;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber
Params.amax = 15; % Initial value for the aperture size
Params.Cs = 0;
Params.df = 0;
%% Sample preparation
ConeR = sqrt(X.^2 + Y.^2); % in Angstrom
BottomR = 2; % in Angstrom
ConeHeight = 500 : 100 : 2500;
SemiAngle = 15 : 30;
MOI = zeros(length(SemiAngle), length(ConeHeight));
% Each column records the MOI data with respect to SemiAngle for a certain
% potential peak value.
figure;
LegendStr = strings([1, length(ConeHeight)]);
for ConeH_Idx = 1 : length(ConeHeight)
    ConePot = -ConeHeight(ConeH_Idx) / BottomR * ConeR + ConeHeight(ConeH_Idx);
    ConePot(ConePot < 0) =0;
    ConeTF = exp(1i * InterCoeff * ConePot / 1e3);
    for SemiAngle_Idx = 1 : length(SemiAngle)
        Params.amax = SemiAngle(SemiAngle_Idx);
        Probe = ProbeCreate(Params, 0, 0, Lx, Ly, Nx, Ny);
        TransWave = Probe .* ConeTF;
        RecWave = ifftshift(fft2(fftshift(TransWave)) * dx * dy);
        DetectInten = abs(RecWave.^2);
        MOI(SemiAngle_Idx, ConeH_Idx) = sum(sum(DetectInten .* (Fx.^2 + Fy.^2) / (Lx * Ly)));
    end
    plot(SemiAngle, MOI( : , ConeH_Idx));
    hold on;
    LegendStr(ConeH_Idx) = [num2str(ConeHeight(ConeH_Idx)), ' V'];
end
legend(LegendStr);
xlabel('Aperture size (mrad)');
ylabel('MOI');