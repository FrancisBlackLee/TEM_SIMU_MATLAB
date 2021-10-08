function SCEM_X(Lx, Ly, params, transFuncs, sliceDists, stackNum, destDir,...
    aberrType, camera, preference, varargin)
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
%   camera -- optional: insert camera or not, syntax: 'camera' or 'none'
%       (default);
%   preference -- 'column'/'row': writing by column (default) or row, note
%       that the aberrType must be specified and the camera should be set 
%       as 'camera' before preference is specified;
%   varargin -- specific parameters to define the writing mode (default:
%       double), note that the preferrence must be specified before
%       varargin are specified;
% Note: scem image data will be saved as *.mat files; camera recordings
%   will be saved as binary files in newly created folders.

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
    destDir = 'tmp_scem_results';
    CreateNewFolder(destDir);
    aberrType = 'reduced';
    camera = 'none';
elseif nargin == 7
    DirectoryExist;
    aberrType = 'reduced';
    camera = 'none';
elseif nargin == 8
    DirectoryExist;
    ValidAberrationType;
    camera = 'none';
elseif nargin == 9
    DirectoryExist;
    ValidAberrationType;
    ValidCamera;
    if strcmp(camera, 'camera')
        preference = 'column';
    end
elseif nargin >= 10
    DirectoryExist;
    ValidAberrationType;
    ValidCamera;
    ValidPreferrence;
end

[Ny, Nx, sliceNum] = size(transFuncs);
wavLen = HighEnergyWavLen_X(params.KeV);
% generate fftshifted Fresnel propagation kernels:
sampleThickness = stackNum * sum(sliceDists);
shiftPropKer = 1i * ones(Ny, Nx, sliceNum);
for sliceIdx = 1 : sliceNum
    shiftPropKer(:, :, sliceIdx) = fftshift(FresnelPropKernel_X(Lx, Ly,...
        Nx, Ny, wavLen, sliceDists(sliceIdx)));
end

% fftshift all the transmission functions in place
transFuncs = fftshift(transFuncs, 1);
transFuncs = fftshift(transFuncs, 2);

dfNum = length(params.dfSeries);
if strcmp(camera, 'none')
    pinholeNum = length(params.pinholeRadii);
    % generate mesh for real-space coordinates for determining pinholes:
    xAxis = InitAxis(Lx, Nx);
    yAxis = InitAxis(Ly, Ny);
    [xMesh, yMesh] = meshgrid(xAxis, yAxis);
end

scanNx = length(params.scanx);
scanNy = length(params.scany);

taskNum = dfNum * scanNy * scanNx;
wbHandle = waitbar(0, 'scanning...');
% Initialize otf
otf = 1i * ones(Ny, Nx);

% fftshift params.lowerAperture for less computation cost
params.lowerAperture = fftshift(params.lowerAperture);

for dfIdx = 1 : dfNum
    df = params.dfSeries(dfIdx);
    dfFolder = sprintf('df=%.4fAngs', df);
    dfFolder = fullfile(destDir, dfFolder);
    CreateNewFolder(dfFolder);
    
    UpdateOTF;
    
    % create a 3D matrix to saved images for various pinhole radii, hope it
    % won't trigger memory overflow:)
    if strcmp(camera, 'none')
        scemImg = zeros(scanNy, scanNx, pinholeNum);
    end
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
            
            secondPropDist = -df - sampleThickness;
            shiftedSecondPropKer = fftshift(FresnelPropKernel_X(Lx, Ly, Nx, Ny,...
                wavLen, secondPropDist));
            
            % real space to reciprocal space:
            wave = fft2(wave);
            wave = shiftedSecondPropKer .* wave;
            wave = params.lowerAperture .* wave;
            % reciprocal space to real space:
            wave = ifft2(wave);
            waveI = abs(wave.^2);
            
            if strcmp(camera, 'none')
                % calculate relative position of the pinhole center to the
                % probe center: relative x
                relativeXMesh = xMesh - params.scanx(xIdx);
                % alter to satisfy periodic boundary condition
                alterIndices = find(relativeXMesh < -Lx / 2.0);
                relativeXMesh(alterIndices) = relativeXMesh(alterIndices) + Lx;

                alterIndices = find(relativeXMesh > Lx / 2.0);
                relativeXMesh(alterIndices) = relativeXMesh(alterIndices) - Lx;

                % relative y
                relativeYMesh = yMesh - params.scany(yIdx);
                % alter to satisfy periodic boundary condition
                alterIndices = find(relativeYMesh < -Ly / 2.0);
                relativeYMesh(alterIndices) = relativeYMesh(alterIndices) + Ly;

                alterIndices = find(relativeYMesh > Ly / 2.0);
                relativeYMesh(alterIndices) = relativeYMesh(alterIndices) - Ly;

                % relative distance
                rMesh = sqrt(relativeXMesh.^2 + relativeYMesh.^2);
                % fftshift rMesh to further reduce computation cost:
                rMesh = fftshift(rMesh);

                for pinholeIdx = 1 : pinholeNum
                    pinhole = (rMesh < params.pinholeRadii(pinholeIdx));
                    scemImg(yIdx, xIdx, pinholeIdx) = sum(waveI .* pinhole, 'all');
                end
            elseif strcmp(camera, 'camera')
                filename = ['intensity_y', num2str(yIdx), '_x', num2str(xIdx), '.bin'];
                filename = fullfile(dfFolder, filename);
                waveI = ifftshift(waveI);
                WriteBinaryFile(filename, waveI, preference, varargin{:});
            end
        end
    end
    
    if strcmp(camera, 'none')
        filename = 'scem_images.mat';
        filename = fullfile(dfFolder, filename);
        save(filename, 'scemImg');
    end
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

    function ValidCamera
        if ~(strcmp(camera, 'none') || strcmp(camera, 'camera'))
            error('Error: invalid input of camera option!');
        end
    end

    function ValidPreferrence
        if ~(strcmp(preference, 'column') || strcmp(preference, 'row'))
            error('Error: invalid input of preferrence!');
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

