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
function [stemImg] = ADF_STEM_X(Lx, Ly, params, TransFuncs, SliceDist, StackNum, CBEDoption, CBEDdir)
%ADF_STEM_X.m is a specially designed multislice interface.
%   Lx, Ly -- sampling sidelength in angstrom;
%   params -- ADF-STEM parameter setting:
%       params.KeV -- beam energy in KeV;
%       params.Cs3 -- 3rd order spherical aberration in mm;
%       params.Cs5 -- 5th order spherical aberration in mm;
%       params.df -- defocus in angstrom (a negative defocus to eliminate
%           the effect of spherical aberrations, different from old versions);
%       params.NA -- numerical aperture in mrad;
%       params.scanx -- coordinate array for x directional scanning;
%       params.scany -- coordinate array for y directional scanning;
%       params.detector -- this detector is an Ny by Nx logical matrix,
%           which does not have to be annular;
%   TransFuncs -- transmission functions stored in a 3D matrix, keep it
%       small to avoid memory overflow;
%   SliceDist -- distance from one slice to the next slice;
%   StackNum -- number of repeations of the input TransFuncs in z
%       direction;
%   CBEDoption -- whether to save the CBED data to an existed empty local
%       folder (default value is 0, when this is not input):
%       CBEDoption = 0: no; CBEDoption = 1: yes;
%   CBEDdir -- an existed empty local folder to save the CBED data.

if(nargin == 6)
    CBEDoption = 0;
else
    if CBEDoption == 1
        if ~isfolder(CBEDdir)
            errorMessage = sprintf('Error: %s does not exist!\n', CBEDdir);
            uiwait(warndlg(errorMessage));
            return;
        end
    end
end
[Ny, Nx, SliceNum] = size(TransFuncs);
dx = Lx / Nx;
dy = Ly / Ny;
WavLen = 12.3986 / sqrt((2 * 511.0 + params.KeV) * params.KeV);
% generate fftshifted Fresnel propagation kernels:
shiftPropKer = zeros(Ny, Nx, SliceNum) + 1i * ones(Ny, Nx, SliceNum);
for SliceIdx = 1 : SliceNum
    shiftPropKer(:, :, SliceIdx) = fftshift(FresnelPropKernel_X(Lx, Ly, Nx, Ny, WavLen, SliceDist(SliceIdx)));
end
% generate objective transfer function with an aperture:
OTF = CircApert_X(Lx, Ly, Nx, Ny, WavLen, params.NA) .* ObjTransFunc_X(params, Lx, Ly, Nx, Ny);

% fftshift all the transmission function in place:
TransFuncs = fftshift(TransFuncs, 1);
TransFuncs = fftshift(TransFuncs, 2);
% start scanning:
scanNx = length(params.scanx);
scanNy = length(params.scany);
stemImg = zeros(scanNy, scanNx);
TotalCompTask = scanNy * scanNx;
process = waitbar(0, 'start scanning');
for iy = 1 : scanNy
    for ix = 1 : scanNx
        DoneRatio = ((iy - 1) * scanNx + ix - 1) / TotalCompTask;
        waitbar(DoneRatio, process, [num2str(roundn(DoneRatio, -3) * 100), '%']);
        TempWave = fftshift(GenerateProbe_X(OTF, params.scanx(ix), params.scany(iy), Lx, Ly, Nx, Ny));
        for StackIdx = 1 : StackNum
            for SliceIdx = 1 : SliceNum
                TempWave = TempWave .* TransFuncs(:,:,SliceIdx);
                TempWave = ifft2(shiftPropKer(:, :, SliceIdx) .* fft2(TempWave));
            end
        end
        TempInten = abs((ifftshift(fft2(TempWave))*dx*dy).^2);
        stemImg(iy, ix) = sum(sum(TempInten .* (params.detector)));
        if CBEDoption == 1
            CBEDname = strcat(CBEDdir, '\cbed_', num2str(iy), '_', num2str(ix), '.txt');
            save(CBEDname, 'TempInten', '-ascii', '-double', '-tabs');
        end
    end
end
delete(process);

end

