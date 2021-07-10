function [limited] = BandwidthLimit(unlimited, Lx, Ly, Nx, Ny, proportion)
%BandwidthLimit limits the input 2D matrix in the frequency space with 
%input cutoff proportion ranging from 0 to 1.
%   unlimited -- matrix to be bandwidth limited;
%   Lx, Ly, Nx, Ny -- sampling parameters;
%   proportion -- bandwidth limit proportion;
%   limited -- bandwidth limited matrix;

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

dx = Lx / Nx;
dy = Ly / Ny;
fx = -1 / (2 * dx) : 1 / Lx : 1 / (2 * dx) - 1 / Lx;
fy = -1 / (2 * dy) : 1 / Ly : 1 / (2 * dy) - 1 / Ly;
[Fx, Fy] = meshgrid(fx, fy);
FreqSquare = Fx.^2 + Fy.^2;
Aperture = ones(size(Fx));
R_apert = proportion * min(1 / (2 * dx), 1 / (2 * dy));
Aperture(FreqSquare > R_apert^2) = 0;
reciprocal = fftshift(Aperture) .* fft2(fftshift(unlimited));
limited = ifftshift(ifft2(reciprocal));
amp = abs(limited);
limited = limited ./ amp;

end

