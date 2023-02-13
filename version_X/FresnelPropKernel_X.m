function [propKernel] = FresnelPropKernel_X(Lx, Ly, Nx, Ny, wavLen, propDist, tiltX, tiltY)
%FresnelPropKernel_X.m computes the Fresnel propagation kernel.
%   Lx, Ly -- sampling side length;
%   Nx, Ny -- sampling number;
%   wavLen -- wavelength;
%   propDist -- propagation distance
%   tiltX -- small tilt angle in mrad in x direction
%   tiltY -- small tilt angle in mrad in y direction
% Note: X denotes an experimental version!

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

fx = InitFreqAxis(Lx, Nx);
fy = InitFreqAxis(Ly, Ny);
[FX, FY] = meshgrid(fx, fy);
if nargin == 6
    propKernel = exp(-1i * pi * wavLen * propDist * (FX.^2 + FY.^2));
elseif nargin == 8
    tx = tiltX * 1.0e-3;
    ty = tiltY * 1.0e-3;
    propKernel = exp(-1i * pi * wavLen * propDist * (FX.^2 + FY.^2) +...
        1i * 2 * pi * propDist * (FX * tan(tx) + FY * tan(ty)));
else
    error('Incorrect number of arguments.');
end

end



