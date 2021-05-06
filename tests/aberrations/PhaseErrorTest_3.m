%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2020  Francis Black Lee and Li Xian

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
% PhaseErrorTest_3.m -- test AberrationPhaseShift_X.m
clc;
clear;
close all;
%   Aberration items:
%       Aberr{1} = [C10, C12];
%       Aberr{2} = [C21, C23];
%       Aberr{3} = [C30, C32, C34];
%       Aberr{4} = [C41, C43, C45];
%       Aberr{5} = [C50, C52, C54, C56];
%       C10  C1  Defocus
%       C12  A1  2-Fold astigmatism
%       C21  B2  Axial coma
%       C23  A2  3-Fold astigmatism
%       C30  C3  3rd order spherical aberration
%       C32  S3  Axial star aberration
%       C34  A3  4-Fold astigmatism
%       C41  B4  4th order axial coma
%       C43  D4  3-Lobe aberration
%       C45  A4  5-Fold astigmatism
%       C50  C5  5th order spherical aberration
%       C52  S5  5th order axial star
%       C54  R5  5th order rosette
%       C56  A5  6-Fold astigmatism
%   Real rotational angle in degree:
%       RealRotAngle{1} = [phi_10, phi_12];
%       RealRotAngle{2} = [phi_21, phi_23];
%       RealRotAngle{3} = [phi_30, phi_32, phi_34];
%       RealRotAngle{4} = [phi_41, phi_43, pni_45];
%       RealRotAngle{5} = [phi_50, phi_52, phi_54, phi_56];
%% Parameter setting:
% unit of aberrations is angstrom
aberration = InitObjectiveLensAberrations_X();
aberration.C1 = 0;
aberration.A1 = 0;
aberration.B2 = 0;
aberration.A2 = 0;
aberration.C3 = 0.057e7;
aberration.S3 = 0;
aberration.A3 = 0;
aberration.A4 = 0;

aberration.A1_angle = 0;
aberration.B2_angle = 0;
aberration.A2_angle = 0;
aberration.S3_angle = 0;
aberration.A3_angle = 0;

dx = 0.05; % Angstrome, expand reciprocal space to about 70 mrad.
dy = dx;
Nx = 2048;
Ny = 2048;
Lx = Nx * dx;
Ly = Ny * dy;
x = -Lx / 2 : dx : Lx / 2 - dx;
y = -Ly / 2 : dy : Ly / 2 - dy;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;

KeV = 200;
wavLen = HighEnergyWavLen_X(KeV);

aperture = CircApert_X(Lx, Ly, Nx, Ny, wavLen, 20);

phaseError = AberrationPhaseShift_X(aberration, wavLen, Lx, Ly, Nx, Ny);
wrappedPhaseError = mat2gray(wrapTo2Pi(phaseError)) .* aperture;

figure;
imshow(wrappedPhaseError, []);
colormap('gray'); axis square;

opticalTransFunc = exp(-1i * phaseError) .* aperture;
probe = GenerateProbe_X(opticalTransFunc, 0, 0, Lx, Ly, Nx, Ny);

probeI = abs(probe.^2);
cropRange = 1024 - 256 + 1 : 1024 + 256;
cropProbeI = probeI(cropRange, cropRange);
figure;
imshow(cropProbeI, []);
colormap('gray'); axis square;

figure;
plot(x(513 : 1536), probeI(1025, 513 : 1536));
xlabel('$x (\AA)$', 'interpreter', 'latex');
ylabel('Intensity');