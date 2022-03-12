function [limited] = BandwidthLimit(unlimited, Lx, Ly, Nx, Ny, proportion, shiftMode)
%BandwidthLimit limits the input 2D matrix in the frequency space with 
%input cutoff proportion ranging from 0 to 1.
%   unlimited -- matrix to be bandwidth limited;
%   Lx, Ly, Nx, Ny -- sampling parameters;
%   proportion -- bandwidth limit proportion;
%   limited -- bandwidth limited matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2022  Francis Black Lee (Li Xian)

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

if nargin == 6
    shiftMode.in = false;
    shiftMode.out = false;
end

dx = Lx / Nx;
dy = Ly / Ny;
fx = InitFreqAxis(Lx, Nx);
fy = InitFreqAxis(Ly, Ny);
[fxMesh, fyMesh] = meshgrid(fx, fy);
freqSquare = fxMesh.^2 + fyMesh.^2;
r = proportion * min(1 / (2 * dx), 1 / (2 * dy));
aperture = (freqSquare <= r^2);

if shiftMode.in
    reciprocal = fftshift(aperture) .* fft2(unlimited);
else
    reciprocal = fftshift(aperture) .* fft2(fftshift(unlimited));
end

if shiftMode.out
    limited = ifft2(reciprocal);
else
    limited = ifftshift(ifft2(reciprocal));
end

amp = abs(limited);
limited = limited ./ amp;

end


