%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADF-STEM sample 0: silicon [110]
clc;
close all;
clear all;
%% Lattice generation: silicon [110]
W_B = waitbar(0, 'Preparing the specimen...');

LattConst = [3.84, 5.43, 0]; % [a b]
LayerDist = [1.9198, 1.9198]; % distance between each slice
CellNum = [6, 4]; % expand the unit cell by Expan_Nx = 3 and Expan_Ny = 2, adaptive
% Laters: Each column for an atom
LayerA = [14, 14; 0, 0.5; 0, 0.75];
LayerB = [14, 14; 0, 0.5; 0.25, 0.5];
%% basic settings
% sampling:
Lx = CellNum(1) * LattConst(1);
Ly = CellNum(2) * LattConst(2);
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
Params.KeV = 200;
InterCoeff = InteractionCoefficient(Params.KeV);
WaveLength = 12.3986 / sqrt((2 * 511.0 + Params.KeV) * Params.KeV);  %wavelength
WaveNumber = 2 * pi / WaveLength;     %wavenumber
Params.amax = 10.37;
Params.Cs = 1.3;
Params.df = 600;

detector = Fx.^2 + Fy.^2;
detector_cri = detector;
HighAngle = 200 * 0.001;
LowAngle = 40 * 0.001;
detector((detector_cri > (sin(LowAngle) / WaveLength)^2) & (detector_cri < (sin(HighAngle) / WaveLength)^2)) = 1;
detector((detector_cri < (sin(LowAngle) / WaveLength)^2) | (detector_cri > (sin(HighAngle) / WaveLength)^2)) = 0;
%% Transmission functions
% Layer A:
Proj_PotA = MultiProjPot_conv_0(LayerA, CellNum, LattConst, Lx, Ly, Nx, Ny);
% Layer B:
Proj_PotB = MultiProjPot_conv_0(LayerB, CellNum, LattConst, Lx, Ly, Nx, Ny);
% test
figure;
imagesc(x, y, Proj_PotA);
colormap('gray');
figure;
imagesc(x, y, Proj_PotB);
colormap('gray');

TF_A = exp(1i * InterCoeff * Proj_PotA / 1000);
TF_B = exp(1i * InterCoeff * Proj_PotB / 1000);
TF_A = BandwidthLimit(TF_A, Lx, Ly, Nx, Ny, 0.67);
TF_B = BandwidthLimit(TF_B, Lx, Ly, Nx, Ny, 0.67);
TransFuncs(:, :, 1) = TF_A;
TransFuncs(:, :, 2) = TF_B;

waitbar(0, W_B, 'Specimen preparation completed, start scanning...');
%% Scanning module
% Scanning parameters:
Scan_Nx = 16; % scanning sampling number, adaptive
Scan_Ny = 16;
Scan_Lx = Lx / 1.5; % scanning side length, adaptive
Scan_Ly = Ly / 1.5;
Scan_dx = Scan_Lx / Scan_Nx;
Scan_dy = Scan_Ly / Scan_Ny;
ADF_x = -Scan_Lx / 2 : Scan_dx : Scan_Lx / 2 - Scan_dx;
ADF_y = -Scan_Ly / 2 : Scan_dy : Scan_Ly / 2 - Scan_dy;
STEM_IMAGE = zeros(Scan_Ny, Scan_Nx);
figure;
StackNum = 20; % determines the thickness of the specimen
TotalNum = Scan_Ny * Scan_Nx;
for i=1:Scan_Ny
    yp = ADF_y(i);
    for j=1:Scan_Nx
        xp = ADF_x(j);
        Probe = ProbeCreate(Params, xp, yp, Lx, Ly, Nx, Ny);
        Trans_Wave = multislice(Probe, WaveLength, Lx, Ly, TransFuncs, LayerDist, StackNum);
        Trans_Wave_Far = ifftshift(fft2(fftshift(Trans_Wave))*dx^2);
        DetectInten = abs(Trans_Wave_Far.^2).*detector;
        STEM_IMAGE(i,j) = sum(sum(DetectInten));
        CurrentNum = (i - 1) * Scan_Nx + j;
        imagesc(ADF_x, ADF_y, STEM_IMAGE);
        map = colormap(gray);
        axis square;
        title('Example');
        drawnow;
        % Save gif to local directory
%         F = getframe(gcf); 
%         I = frame2im(F); 
%         [I, map] = rgb2ind(I, 256); 
%         if CurrentNum == 1 
%             imwrite(I,map,'D:\Francis. B. Lee\Practice\Conventional Multislice in MATLAB\Specimen_Thickness\Secret\secret.gif','gif','Loopcount',inf,'DelayTime',0.02); 
%         else
%             imwrite(I,map,'D:\Francis. B. Lee\Practice\Conventional Multislice in MATLAB\Specimen_Thickness\Secret\secret.gif','gif','WriteMode','append','DelayTime',0.02); 
%         end
        waitbar(roundn(CurrentNum / TotalNum, -3), W_B, [num2str(roundn((CurrentNum / TotalNum), -3) * 100), '%']);
    end
end
delete(W_B);