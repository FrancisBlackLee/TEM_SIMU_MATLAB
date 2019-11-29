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
function [ProjPot] = MonoProjPot_conv_0(AtomType, ScaleCoord, CellNum, LattConst, Lx, Ly, Nx, Ny)
%MonoProjPot_conv_0.m calculates the projected potential for one type of
%atom, more comprehensive function to be released soon.
%   ScaleCoord -- scaled planar coordinates, syntax: [ScaleX1,..., ScaleXN;
%       ScaleY1, ScaleYN];
%   CellNum -- number of unit cells to be included, sytax: [CellNumX, CellNumY];
%   LattConst -- lattice constants, syntax: [a, b];
%   Lx, Ly, Nx, Ny -- sampling parameters;

AtomNum = size(ScaleCoord, 2);
SingPot_fft = fft2(fftshift(ProjectedPotential(Lx, Ly, Nx, Ny, AtomType, 0, 0)));
% Sampling:
dx = Lx / Nx;
dy = Ly / Ny;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
% Build the convolution kernel:
Kernel = 0;
ScaleXshift = CellNum(1) / 2;
ScaleYshift = CellNum(2) / 2;
for Atom_Idx = 1 : AtomNum
    for Cellx_Idx = 0 : CellNum(1) - 1
        CoordX = LattConst(1) * (ScaleCoord(1, Atom_Idx) + Cellx_Idx - ScaleXshift);
        for Celly_Idx = 0 : CellNum(2) - 1
            CoordY = LattConst(2) * (ScaleCoord(2, Atom_Idx) + Celly_Idx - ScaleYshift);
            Kernel = Kernel + exp(-1i * 2 * pi *(Fx * CoordX + Fy * CoordY));
        end
    end
end
Kernel = fftshift(Kernel);
ProjPot = real(ifftshift(ifft2(Kernel .* SingPot_fft)));

end

