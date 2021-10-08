function SCEM_Preprocessing_X(Lx, Ly, params, transFuncs, sliceDists,...
    stackNum, destDir, aberrType, preference, varargin)
%SCEM_Preprocessing_X.m simulates the beam-specimen interaction using
%multislice method under scanning confocal mode, and outputs the wave
%function at each scanning point to the destination directory.
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
%   preferrence -- 'column'/'row': writing by column (default) or row, note
%       that the aberrType must be specified before preferrence is
%       specified;
%   varargin -- specific parameters to define the writing mode (default:
%       double), note that the preferrence must be specified before
%       varargin are specified;

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
    preference = 'column';
elseif nargin == 7
    DirectoryExist;
    aberrType = 'reduced';
    preference = 'column';
elseif nargin == 8
    DirectoryExist;
    ValidAberrationType;
    preference = 'column';
elseif nargin >= 9
    DirectoryExist;
    ValidAberrationType;
    ValidPreferrence;
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
scanNx = length(params.scanx);
scanNy = length(params.scany);

taskNum = dfNum * scanNy;
wbHandle = waitbar(0, 'scanning...');
% Initialize otf
otf = 1i * ones(Ny, Nx);
for dfIdx = 1 : dfNum
    df = params.dfSeries(dfIdx);
    dfFolder = ['df=', num2str(df, '%.4f'), 'Angs'];
    dfFolder = fullfile(destDir, dfFolder);
    if(~CreateNewFolder(dfFolder))
        return;
    end
    
    UpdateOTF;
    
    for yIdx = 1 : scanNy
        doneRatio = ((dfIdx - 1) * scanNy + yIdx - 1) / taskNum;
        wbMessage = sprintf('df: %d / %d, line: %d / %d completed',...
            dfIdx - 1, dfNum, yIdx - 1, scanNy);
        waitbar(doneRatio, wbHandle, wbMessage);
        for xIdx = 1 : scanNx
            wave = fftshift(GenerateProbe_X(otf, params.scanx(xIdx),...
                params.scany(yIdx), Lx, Ly, Nx, Ny));
            for stackIdx = 1 : stackNum
                for sliceIdx = 1 : sliceNum
                    wave = wave .* transFuncs(:, :, sliceIdx);
                    wave = ifft2(shiftPropKer(:, :, sliceIdx) .* fft2(wave));
                end
            end
            
            wave = ifftshift(wave);
            % save wave function to dfFolder
            filename = ['wave_y', num2str(yIdx), '_x', num2str(xIdx), '.bin'];
            filename = fullfile(dfFolder, filename);
            WriteComplexBinaryFile(filename, wave, preference, varargin{:});
        end
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

    function ValidPreferrence
        if ~(strcmp(preference, 'column') || strcmp(preference, 'row'))
            error('Error: invalid input of preferrence!');
        end
    end
    
    function UpdateOTF
        if strcmp(aberrType, 'reduced')
            params.df = df;
            otf = params.aperture .* ObjTransFunc_X(params, Lx, Ly, Nx, Ny);
        elseif strcmp(aberrType, 'full')
            params.aberration.C1 = df;
            otfPhase = AberrationPhaseShift_X(params.aberration,...
                wavLen, Lx, Ly, Nx, Ny);
            otf = params.aperture .* exp(-1i * otfPhase);
        end
    end

end

