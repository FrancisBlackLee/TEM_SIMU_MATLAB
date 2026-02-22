function [stemImg, seImg] = STEM_X_SE_var_z_gpu(Lx, Ly, params, transFuncs, seObjFuncs,...
    mfp, sliceDist, stackNum, z0, aberrType)
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
%   seObjFuncs -- secondary electron object functions, 3D matrix like
%       transFuncs;
%   mfp -- mean free path
%   sliceDist -- distance from one slice to the next slice;
%   stackNum -- number of repeations of the input TransFuncs in z
%       direction;
%   aberrType -- aberration type: 'reduced' (C3, C5, df) or 'full' (up to
%       5th order aberrations and their real-space angles (must be
%       specified when cbedOption is specified, otherwise the function will
%       exit with an error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Xinyuan Zhang, Francis Black Lee and Li Xian

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

if nargin == 9
    aberrType = 'reduced';
end
[Ny, Nx, sliceNum] = size(transFuncs);
dx = Lx / Nx;
dy = Ly / Ny;
wavLen = HighEnergyWavLen_X(params.KeV);
% generate fftshifted Fresnel propagation kernels:
shiftPropKer = 1i * ones(Ny, Nx, sliceNum, "single", "gpuArray");
for sliceIdx = 1 : sliceNum
    shiftPropKer(:, :, sliceIdx) = fftshift(gpuArray(single(FresnelPropKernel_X(Lx, Ly,...
        Nx, Ny, wavLen, sliceDist(sliceIdx)))));
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

% fftshift all the secondary electron object functions in place:
seObjFuncs = fftshift(seObjFuncs, 1);
seObjFuncs = fftshift(seObjFuncs, 2);

% start scanning:
scanNx = length(params.scanx);
scanNy = length(params.scany);
detectorNum = size(params.detector, 3);
stemImg = zeros(scanNy, scanNx, detectorNum, "single", "gpuArray");
seImg = zeros(scanNy, scanNx, "single", "gpuArray");

process = waitbar(0, 'start scanning');
for iy = 1 : scanNy
    for ix = 1 : scanNx
        tempWave = fftshift(GenerateProbe_X(otf, params.scanx(ix), params.scany(iy),...
            Lx, Ly, Nx, Ny));
        
        depth = 0.0;
        seSignal = 0.0;
        for stackIdx = 1 : stackNum
            for sliceIdx = 1 : sliceNum
                tempWave = tempWave .* transFuncs(:,:,sliceIdx);
                tempWave = ifft2(shiftPropKer(:, :, sliceIdx) .* fft2(tempWave));
                depth = depth + sliceDist(sliceIdx);
                
                generationRate = sum((abs(tempWave.^2) .* seObjFuncs(:, :, sliceIdx)), 'all');
                seTransport = exp(-(depth - z0(iy, ix)) / mfp(iy, ix));
                seSignal = seSignal + generationRate * seTransport;
            end
        end
        tmpCbed = abs((ifftshift(fft2(tempWave))).^2);
        
        for detectorIdx = 1 : detectorNum
            stemImg(iy, ix, detectorIdx) = sum(tmpCbed .*...
                params.detector(:, :, detectorIdx), 'all');
        end
        
        seImg(iy, ix) = seSignal;
    end

    doneRatio = iy / scanNy;
    waitbar(doneRatio, process, [num2str(roundn(doneRatio, -3) * 100), '%']);
end
delete(process);

end

