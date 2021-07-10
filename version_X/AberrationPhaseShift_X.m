function [phaseShift] = AberrationPhaseShift_X(aberration, wavLen, Lx, Ly, Nx, Ny)
%AberrationPhaseShift_X() converts aberration items to phase shift.
% Input:
%   Aberration items:
%       C1  Defocus
%       A1  2-Fold astigmatism
%       B2  Axial coma
%       A2  3-Fold astigmatism
%       C3  3rd order spherical aberration
%       S3  Axial star aberration
%       A3  4-Fold astigmatism
%       B4  4th order axial coma
%       D4  3-Lobe aberration
%       A4  5-Fold astigmatism
%       C5  5th order spherical aberration
%       S5  5th order axial star
%       R5  5th order rosette
%       A5  6-Fold astigmatism
%
%       A1_angle  2-Fold astigmatism
%       B2_angle  Axial coma
%       A2_angle  3-Fold astigmatism
%       S3_angle  Axial star aberration
%       A3_angle  4-Fold astigmatism
%       B4_angle  4th order axial coma
%       D4_angle  3-Lobe aberration
%       A4_angle  5-Fold astigmatism
%       S5_angle  5th order axial star
%       R5_angle  5th order rosette
%       A5_angle  6-Fold astigmatism
%   Lx, Ly, Nx, Ny -- sampling parameters;
%   Version note: X denotes an experimental version. For more details,
%       please refer to Sec. 2.8 of Advanced Computig Electron Microscopy.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

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

aberr{1} = [aberration.C1, aberration.A1];
aberr{2} = [aberration.B2, aberration.A2];
aberr{3} = [aberration.C3, aberration.S3, aberration.A3];
aberr{4} = [aberration.B4, aberration.D4, aberration.A4];
aberr{5} = [aberration.C5, aberration.S5, aberration.R5, aberration.A5];

realRotAngle{1} = [0, aberration.A1_angle];
realRotAngle{2} = [aberration.B2_angle, aberration.A2_angle];
realRotAngle{3} = [0, aberration.S3_angle, aberration.A3_angle];
realRotAngle{4} = [aberration.B4_angle, aberration.D4_angle, aberration.A4_angle];
realRotAngle{5} = [0, aberration.S5_angle, aberration.R5_angle, aberration.A5_angle];

phaseShift = MultiAberrPhaseError_X_v2(aberr, realRotAngle, wavLen, Lx, Ly, Nx, Ny);

end

