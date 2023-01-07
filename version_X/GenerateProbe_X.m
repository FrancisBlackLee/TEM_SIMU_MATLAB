function [probe] = GenerateProbe_X(varargin)
%GenerateProbe_X.m generates an electron probe.
% GenerateProbe_X(otf, xp, yp, Lx, Ly, Nx, Ny):
%   OTF -- prepared objective transfer function in reciprocal space, note
%       that ObjTransFunc_X.m does not generate the OTF with an aperture,
%       thus the input OTF must have been multiplied by an aperture in
%       advance;
%   Lx, Ly, Nx, Ny -- sampling parameters, L denotes side length and N the
%       sampling number in real space;
%   xp, yp -- probe position in real space;
%
% GenerateProbe_X(params):
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

if nargin == 7
    otf = varargin{1};
    xp = varargin{2};
    yp = varargin{3};
    Lx = varargin{4};
    Ly = varargin{5};
    Nx = varargin{6};
    Ny = varargin{7};

elseif nargin == 1
    params = varargin{1};
    xp = params.xp;
    yp = params.yp;
    Lx = params.Lx;
    Ly = params.Ly;
    Nx = params.Nx;
    Ny = params.Ny;
    wavLen = HighEnergyWavLen_X(params.KeV);
    if strcmp(params.type, 'reduced')
        otf = params.aperture .* ObjTransFunc_X(params, Lx, Ly, Nx, Ny);
    elseif strcmp(params.type, 'full')
        otfPhase = AberrationPhaseShift_X(params.aberration, wavLen,...
            Lx, Ly, Nx, Ny);
        otf = params.aperture .* exp(-1i * otfPhase);
    else
        error('Invalid aberration type!\n');
    end
else
    error("Incorrect number of input arguments.");
end

dx = Lx / Nx;
dy = Ly / Ny;
fx = InitFreqAxis(Lx, Nx);
fy = InitFreqAxis(Ly, Ny);
[FX, FY] = meshgrid(fx, fy);
probe = ifftshift(ifft2(fftshift(otf .* exp(-1i * 2 * pi * (FX * xp + FY * yp)))));

normCoeff = sqrt(sum(abs(probe.^2), 'all') * dx * dy);
probe = probe / normCoeff;

end



