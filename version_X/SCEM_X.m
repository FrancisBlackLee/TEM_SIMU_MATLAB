function SCEM_X(Lx, Ly, params, transFuncs, sliceDists, stackNum, destDir,...
    aberrType)
%SCEM_X.m simulates scanning confocal electron microscopy, in which the
%beam-specimen interaction is simulated using multislice method.
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
%       params.upperAperture -- the aperture at the position of the upper
%           objective lens;
%       params.lowerAperture -- the aperture at the position of the lower
%           objective lens;
%       params.pinholeRadii -- the radii of real-space pinholes, which are
%           positioned below the lower objective aperture, the physical
%           unit of radii should coincide with that of Lx and Ly;
%       params.scanx -- coordinate array for x directional scanning;
%       params.scany -- coordinate array for y directional scanning;
%       params.dfSeries -- objective defocus range for SCEM imaging (in
%           Angstrom);
%   transFuncs -- transmission functions stored in a 3D matrix, keep it
%       small to avoid memory overflow;
%   sliceDists -- distance from one slice to the next slice;
%   stackNum -- number of repeations of the input TransFuncs in z
%       direction;
%   destDir -- a local directory to save the exit wave beneath the specimen
%       while scanning the probe with various defocus;
%   aberrType -- aberration type: 'reduced' (C3, C5, df) or 'full' (up to
%       5th order aberrations and their real-space angles;
% Note: data will be saved as *.mat files.

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
    % create a new folder to save the results when destDir is not specified
    destDir = 'tmp_scem_wave';
    CreateNewFolder(destDir);
    aberrType = 'reduced';
elseif nargin == 7
    DirectoryExist;
    aberrType = 'reduced';
elseif nargin == 8
    DirectoryExist;
    ValidAberrationType;
end

[Ny, Nx, sliceNum] = size(transFuncs);
wavLen = HighEnergyWavLen_X(params.KeV);
% generate fftshifted Fresnel propagation kernels:
shiftPropKer = 1i * ones(Ny, Nx, sliceNum);
for sliceIdx = 1 : sliceNum
    shiftPropKer(:, :, sliceIdx) = fftshift(FresnelPropKernel_X(Lx, Ly,...
        Nx, Ny, wavLen, sliceDists(sliceIdx)));
end

% fftshift all the transmission functions in place
transFuncs = fftshift(transFuncs, 1);
transFuncs = fftshift(transFuncs, 2);

dfNum = length(params.dfSeries);
pinholeNum = length(params.pinholeRadii);
scanNx = length(params.scanx);
scanNy = length(params.scany);

taskNum = dfNum * scanNy * scanNx;
wbHandle = waitbar(0, 'scanning...');
% Initialize otf
otf = 1i * ones(Ny, Nx);

% fftshift params.lowerAperture for less computation cost
params.lowerAperture = fftshift(params.lowerAperture);

% generate mesh for real-space coordinates for determining pinholes:
xAxis = InitAxis(Lx, Nx);
yAxis = InitAxis(Ly, Ny);
[xMesh, yMesh] = meshgrid(xAxis, yAxis);
rMesh = sqrt(xMesh.^2 + yMesh.^2);
% fftshift rMesh to further reduce computation cost:
rMesh = fftshift(rMesh);

for dfIdx = 1 : dfNum
    df = params.dfSeries(dfIdx);
    dfFolder = sprintf('df=%.4fAngs', df);
    dfFolder = fullfile(destDir, dfFolder);
    CreateNewFolder(dfFolder);
    
    UpdateOTF;
    
    % create a 3D matrix to saved images for various pinhole radii, hope it
    % won't trigger memory overflow:)
    scemImg = zeros(scanNy, scanNx, pinholeNum);
    for yIdx = 1 : scanNy
        for xIdx = 1 : scanNx
            % update process
            doneRatio = ((dfIdx - 1) * scanNy * scanNx +...
                (yIdx - 1) * scanNx + xIdx - 1) / taskNum;
            wbMessage = sprintf('process %.1f%%, df: %d / %d, line: %d / %d completed',...
                100 * doneRatio, dfIdx - 1, dfNum, yIdx - 1, scanNy);
            waitbar(doneRatio, wbHandle, wbMessage);
            
            wave = fftshift(GenerateProbe_X(otf, params.scanx(xIdx),...
                params.scany(yIdx), Lx, Ly, Nx, Ny));
            for stackIdx = 1 : stackNum
                for sliceIdx = 1 : sliceNum
                    wave = wave .* transFuncs(:, :, sliceIdx);
                    wave = ifft2(shiftPropKer(:, :, sliceIdx) .* fft2(wave));
                end
            end
            
            % real space to reciprocal space:
            wave = fft2(wave);
            wave = params.lowerAperture .* wave;
            % reciprocal space to real space:
            wave = ifft2(wave);
            waveI = abs(wave.^2);
            
            for pinholeIdx = 1 : pinholeNum
                pinhole = (rMesh < params.pinholeRadii(pinholeIdx));
                scemImg(yIdx, xIdx, pinholeIdx) = sum(waveI .* pinhole, 'all');
            end
        end
    end
    
    filename = 'scem_images.mat';
    filename = fullfile(dfFolder, filename);
    save(filename, 'scemImg');
end

delete(wbHandle);

% nested functions:
    function CreateNewFolder(folderName)
        status = mkdir(folderName);
        if status == 0
            error('Error: creating %s failed', folderName);
        end
    end

    function DirectoryExist
        if ~isfolder(destDir)
            error('Error: %s does not exist!\n', destDir);
        end
    end

    function ValidAberrationType
        if ~(strcmp(aberrType, 'reduced') || strcmp(aberrType, 'full'))
            error('Error: invalid input of aberration type!');
        end
    end
    
    function UpdateOTF
        if strcmp(aberrType, 'reduced')
            params.df = df;
            otf = params.upperAperture .* ObjTransFunc_X(params, Lx, Ly, Nx, Ny);
        elseif strcmp(aberrType, 'full')
            params.aberration.C1 = df;
            otfPhase = AberrationPhaseShift_X(params.aberration,...
                wavLen, Lx, Ly, Nx, Ny);
            otf = params.upperAperture .* exp(-1i * otfPhase);
        end
    end

end

