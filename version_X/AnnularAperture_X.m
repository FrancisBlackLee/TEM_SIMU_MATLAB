function [aperture] = AnnularAperture_X(Lx, Ly, Nx, Ny, wavLen,...
    innerAngle, outerAngle)
%AnnularAperture_X.m generates an annular aperture in reciprocal space.
%   Lx, Ly, Nx, Ny -- real-space sampling parameters;
%   wavLen -- wavelength;
%   innerAngle, outerAngle -- inner and outer convergent angles of the
%       aperture (in mrad);
% Note: X denotes an experimental version!

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

fx = InitFreqAxis(Lx, Nx);
fy = InitFreqAxis(Ly, Ny);
[FX, FY] = meshgrid(fx, fy);
angleMesh = sqrt((FX.^2 + FY.^2) * wavLen^2) * 1.0e3;

aperture = (angleMesh > innerAngle) & (angleMesh < outerAngle);

end

