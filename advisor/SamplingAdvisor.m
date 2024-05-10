function [r] = SamplingAdvisor(keV, Lx, Ly, Nx, Ny, bwl)
%SamplingAdvisor.m calculates the pixel size in real space
%and pixel size and maximum angle in reciprocal space.
% Input:
%   keV -- electron kinetic energy in keV;
%   Lx, Ly -- side lengths;
%   Nx, Ny -- pixel numbers;
%   bwl -- bandwidth limit fractor, which is usually 2/3.
% Output:
%   r -- structure containing all the results:
%   r.dx, r.dy -- pixel size in real space in Angs.;
%   r.dax, r.day -- pixel size in reciprocal space in mrad;
%   r.max, r.may -- maximum angles in reciprocal space in mrad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2023  Francis Black Lee (Li Xian)

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

wavLen = HighEnergyWavLen_X(keV);
r.dx = Lx / Nx;
r.dy = Ly / Ny;
r.dax = wavLen * 1e3 / Lx;
r.day = wavLen * 1e3 / Ly;
r.max = wavLen * 1e3 * bwl / (2 * r.dx);
r.may = wavLen * 1e3 * bwl / (2 * r.dy);

end