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
function [PhaseError] = MultiAberrPhaseError_X_v2(Aberr, RealRotAngle,...
    WavLen, Lx, Ly, Nx, Ny)
%MultiAberrPhaseError_X computes the phase error, namely the chi (Greek letter)
%with multiple aberrations.
%   Requests of physical units: all the units for length is Angstrom, and
%   units of aberrations are Angstrom as well.
%   Aberr -- aberration list:
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
%   RealRotAngle -- real rotational angle in degree:
%       RealRotAngle{1} = [phi_10, phi_12];
%       RealRotAngle{2} = [phi_21, phi_23];
%       RealRotAngle{3} = [phi_30, phi_32, phi_34];
%       RealRotAngle{4} = [phi_41, phi_43, pni_45];
%       RealRotAngle{5} = [phi_50, phi_52, phi_54, phi_56];
%   WavLen -- electron wavelength;
%   Lx, Ly, Nx, Ny -- sampling parameters.
%   Version note: X denotes an experimental version. For more details,
%       please refer to Sec. 2.8 of Advanced Computig Electron Microscopy.

AberrNum = [2, 2, 3, 3, 4]; % number of aberrations in each aberration array

dx = Lx / Nx;
dy = Ly / Ny;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[FX, FY] = meshgrid(fx, fy);
PolAng = WavLen * sqrt(FX.^2 + FY.^2);
AziAng = atan2(FY, FX);
% Aberr{n}(Idx) = Cnm, m = (Idx - 1) * 2 + (n - 1) % 2
PhaseError = 0;
for n = 1 : 5
    for Idx = 1 : AberrNum(n)
        m = (Idx - 1) * 2 + mod(n - 1, 2);
        PhaseError = PhaseError + Aberr{n}(Idx) * PolAng.^(n + 1)...
            .* cos(m * (AziAng - deg2rad(RealRotAngle{n}(Idx)))) / (n + 1);
    end
end
PhaseError = 2 * pi / WavLen * PhaseError;

end

