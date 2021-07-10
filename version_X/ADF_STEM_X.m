function [stemImg] = ADF_STEM_X(Lx, Ly, params, transFuncs, sliceDist,...
    stackNum, cbedOption, cbedDir)
%ADF_STEM_X.m is a specially designed multislice interface.
%   Lx, Ly -- sampling sidelength in angstrom;
%   params -- ADF-STEM parameter setting:
%       params.KeV -- beam energy in KeV;
%       params.Cs3 -- 3rd order spherical aberration in mm;
%       params.Cs5 -- 5th order spherical aberration in mm;
%       params.df -- defocus in angstrom (a negative defocus to eliminate
%           the effect of spherical aberrations, different from old versions);
%       params.aperture -- numerical aperture;
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

if(nargin == 6)
    cbedOption = 0;
else
    if cbedOption == 1
        if ~isfolder(cbedDir)
            errorMessage = sprintf('Error: %s does not exist!\n', cbedDir);
            uiwait(warndlg(errorMessage));
            return;
        end
    end
end
[Ny, Nx, sliceNum] = size(transFuncs);
dx = Lx / Nx;
dy = Ly / Ny;
wavLen = 12.3986 / sqrt((2 * 511.0 + params.KeV) * params.KeV);
% generate fftshifted Fresnel propagation kernels:
shiftPropKer = zeros(Ny, Nx, sliceNum) + 1i * ones(Ny, Nx, sliceNum);
for sliceIdx = 1 : sliceNum
    shiftPropKer(:, :, sliceIdx) = fftshift(FresnelPropKernel_X(Lx, Ly,...
        Nx, Ny, wavLen, sliceDist(sliceIdx)));
end
% generate objective transfer function with an aperture:
otf = params.aperture .* ObjTransFunc_X(params, Lx, Ly, Nx, Ny);

% fftshift all the transmission function in place:
transFuncs = fftshift(transFuncs, 1);
transFuncs = fftshift(transFuncs, 2);
% start scanning:
scanNx = length(params.scanx);
scanNy = length(params.scany);
stemImg = zeros(scanNy, scanNx);
totalCompTask = scanNy * scanNx;
process = waitbar(0, 'start scanning');
for iy = 1 : scanNy
    for ix = 1 : scanNx
        doneRatio = ((iy - 1) * scanNx + ix - 1) / totalCompTask;
        waitbar(doneRatio, process, [num2str(roundn(doneRatio, -3) * 100), '%']);
        tempWave = fftshift(GenerateProbe_X(otf, params.scanx(ix), params.scany(iy), Lx, Ly, Nx, Ny));
        for stackIdx = 1 : stackNum
            for sliceIdx = 1 : sliceNum
                tempWave = tempWave .* transFuncs(:,:,sliceIdx);
                tempWave = ifft2(shiftPropKer(:, :, sliceIdx) .* fft2(tempWave));
            end
        end
        tempInten = abs((ifftshift(fft2(tempWave))*dx*dy).^2);
        stemImg(iy, ix) = sum(sum(tempInten .* (params.detector)));
        if cbedOption == 1
            cbedName = strcat(cbedDir, '\cbed_', num2str(iy), '_', num2str(ix), '.txt');
            save(cbedName, 'tempInten', '-ascii', '-double', '-tabs');
        end
    end
end
delete(process);

end

