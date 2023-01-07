function [params] = InitProbeParams_X()
%InitProbeParams_X.m initializes params that GenerateProbe_X(params)
%utilizes.
% Output:
%   params -- parameters for generating the probe:
%       params.KeV -- beam energy in KeV;
%       params.type -- aberration type: 'reduced' (C3, C5, df) or 'full' 
%           (up to 5th order aberrations and their real-space angles);
%       params.Cs3 (optional) -- 3rd order spherical aberration in mm;
%       params.Cs5 (optional) -- 5th order spherical aberration in mm;
%       params.df (optional) -- defocus in angstrom (a negative defocus to 
%           eliminate the effect of spherical aberrations, different from 
%           old versions);
%       params.aberration (optional) -- a more complete set of aberration
%           up to 5th order, it should be initialized with 
%           InitObjectiveLensAberrations_X() and then modified;
%       params.aperture -- numerical aperture;
%       params.xp -- coordinate x of the probe;
%       params.yp -- coordinate y of the probe;
%       params.Lx -- side length for dimension x;
%       params.Ly -- side length for dimension y;
%       params.Nx -- pixel number for dimension x;
%       params.Ny -- pixel number for dimension y;
%
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

params.KeV = 300;
params.type = 'reduced';
params.Cs3 = 0;
params.Cs5 = 0;
params.df = 0;
params.aberration = InitObjectiveLensAberrations_X();
params.aperture = 0;
params.xp = 0;
params.yp = 0;
params.Lx = 10;
params.Ly = 10;
params.Nx = 512;
params.Ny = 512;

end
