function [stemImg] = STEM_X(Lx, Ly, params, transFuncs, sliceDist,...
    stackNum, aberrType, cbedOption, cbedDir, preferrence, varargin)
%STEM_X.m is a specially designed multislice interface for STEM simulation.
%   Lx, Ly -- sampling sidelength in angstrom;
%   params -- STEM parameter setting:
%       params.KeV -- beam energy in KeV;
%       params.Cs3 (optional) -- 3rd order spherical aberration in mm;
%       params.Cs5 (optional) -- 5th order spherical aberration in mm;
%       params.df (optional) -- defocus in angstrom (a negative defocus to 
%           eliminate the effect of spherical aberrations, different from 
%           old versions);
%       params.aberration (optional) -- a more complete set of aberration
%           up to 5th order, it should be initialized with 
%           InitObjectiveLensAberrations_X() and then modified;
%       params.aperture -- numerical aperture;
%       params.scanx -- coordinate array for x directional scanning;
%       params.scany -- coordinate array for y directional scanning;
%       params.detector -- this detector is an Ny by Nx by detectorNum 
%           logical matrix, which do not have to be annular;
%   transFuncs -- transmission functions stored in a 3D matrix, keep it
%       small to avoid memory overflow;
%   sliceDist -- distance from one slice to the next slice;
%   stackNum -- number of repeations of the input TransFuncs in z
%       direction;
%   aberrType -- aberration type: 'reduced' (C3, C5, df) or 'full' (up to
%       5th order aberrations and their real-space angles (must be
%       specified when cbedOption is specified, otherwise the function will
%       exit with an error);
%   cbedOption -- whether to save the CBED data to an existed empty local
%       folder (default value is 0, when this is not input):
%       CBEDoption = 0: no; CBEDoption = 1: yes;
%   cbedDir -- an existed empty local folder to save the CBED data.
%   preferrence -- 'column'/'row': writing by column (default) or row;
%   varargin -- specific parameters to define the writing mode (default:
%       double);

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

if nargin == 6
    cbedOption = 0;
    aberrType = 'reduced';
    preferrence = 'column';
elseif nargin == 7
    cbedOption = 0;
    preferrence = 'column';
elseif nargin == 8
    % when cbedOption == 1 and no input cbedDir, create a temporary one
    if cbedOption == 1
        cbedDir = 'tmp_cbed';
        status = mkdir(cbedDir);
        if status == 0
            error('Error: creating cbed directory failed!');
        end
    end
    preferrence = 'column';
elseif nargin == 9
    if (cbedOption == 1) && (~(strcmp(cbedDir, 'column') || strcmp(cbedDir, 'row')))
        if ~isfolder(cbedDir)
            error('Error: %s does not exist!\n', cbedDir);
        end
    elseif  (cbedOption == 1) && (strcmp(cbedDir, 'column') || strcmp(cbedDir, 'row'))
        % output the cbed patterns without specifying the output directory
        % and with specified binary writing preferrence, create a temporary
        % folder as the output directory and assign the value of the 9th
        % argument to preferrence
        preferrence = cbedDir;
        cbedDir = 'tmp_cbed';
        status = mkdir(cbedDir);
        if status == 0
            error('Error: creating cbed directory failed!');
        end
    end
end
[Ny, Nx, sliceNum] = size(transFuncs);
dx = Lx / Nx;
dy = Ly / Ny;
wavLen = HighEnergyWavLen_X(params.KeV);
% generate fftshifted Fresnel propagation kernels:
shiftPropKer = 1i * ones(Ny, Nx, sliceNum);
for sliceIdx = 1 : sliceNum
    shiftPropKer(:, :, sliceIdx) = fftshift(FresnelPropKernel_X(Lx, Ly,...
        Nx, Ny, wavLen, sliceDist(sliceIdx)));
end
% generate objective transfer function with an aperture:
if strcmp(aberrType, 'reduced')
    otf = params.aperture .* ObjTransFunc_X(params, Lx, Ly, Nx, Ny);
elseif strcmp(aberrType, 'full')
    otfPhase = AberrationPhaseShift_X(params.aberration, wavLen, Lx, Ly, Nx, Ny);
    otf = params.aperture .* exp(-1i * otfPhase);
else
    error('Invalid aberration type!\n');
end

% fftshift all the transmission function in place:
transFuncs = fftshift(transFuncs, 1);
transFuncs = fftshift(transFuncs, 2);
% start scanning:
scanNx = length(params.scanx);
scanNy = length(params.scany);
detectorNum = size(params.detector, 3);
stemImg = zeros(scanNy, scanNx, detectorNum);

totalCompTask = scanNy * scanNx;
process = waitbar(0, 'start scanning');
for iy = 1 : scanNy
    for ix = 1 : scanNx
        doneRatio = ((iy - 1) * scanNx + ix - 1) / totalCompTask;
        waitbar(doneRatio, process, [num2str(roundn(doneRatio, -3) * 100), '%']);
        tempWave = fftshift(GenerateProbe_X(otf, params.scanx(ix), params.scany(iy),...
            Lx, Ly, Nx, Ny));
        for stackIdx = 1 : stackNum
            for sliceIdx = 1 : sliceNum
                tempWave = tempWave .* transFuncs(:,:,sliceIdx);
                tempWave = ifft2(shiftPropKer(:, :, sliceIdx) .* fft2(tempWave));
            end
        end
        tmpCbed = abs((ifftshift(fft2(tempWave))*dx*dy).^2);
        
        for detectorIdx = 1 : detectorNum
            stemImg(iy, ix, detectorIdx) = sum(tmpCbed .*...
                params.detector(:, :, detectorIdx), 'all');
        end
        
        if cbedOption == 1
            cbedName = ['cbed_y', num2str(iy), '_x', num2str(ix), '.bin'];
            cbedName = fullfile(cbedDir, cbedName);
            WriteBinaryFile(cbedName, tmpCbed, preferrence, varargin{:});
        end
    end
end
delete(process);

end

