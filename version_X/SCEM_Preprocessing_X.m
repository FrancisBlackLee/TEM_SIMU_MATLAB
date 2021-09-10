function SCEM_Preprocessing_X(Lx, Ly, params, transFuncs, sliceDist,...
    stackNum, destDir, aberrType, preferrence, varargin)
%UNTITLED Summary of this function goes here
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
%       params.dfSeries -- objective defocus range for SCEM imaging (in
%           Angstrom);
%       params.scany -- coordinate array for y directional scanning;
%   transFuncs -- transmission functions stored in a 3D matrix, keep it
%       small to avoid memory overflow;
%   sliceDist -- distance from one slice to the next slice;
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
    destDir = 'tmp_scem_wave';
    if(~CreateNewFolder(destDir))
        return;
    end
    aberrType = 'reduced';
    preferrence = 'column';
elseif nargin == 7
    if(~DirectoryExist)
        return;
    end
    aberrType = 'reduced';
    preferrence = 'column';
elseif nargin == 8
    if(~DirectoryExist)
        return;
    end
    if(~ValidAberrationType)
        return;
    end
    preferrence = 'column';
elseif nargin == 9
    if(~DirectoryExist)
        return;
    end
    if(~ValidAberrationType)
        return;
    end
    if(~ValidPreferrence)
        return;
    end
end

[Ny, Nx, sliceNum] = size(transFuncs);
wavLen = HighEnergyWavLen_X(params.KeV);
% generate fftshifted Fresnel propagation kernels:
shiftPropKer = 1i * ones(Ny, Nx, sliceNum);
for sliceIdx = 1 : sliceNum
    shiftPropKer(:, :, sliceIdx) = fftshift(FresnelPropKernel_X(Lx, Ly,...
        Nx, Ny, wavLen, sliceDist(sliceIdx)));
end


dfNum = length(params.dfSeries);
scanNx = length(params.scanx);
scanNy = length(params.scany);

taskNum = dfNum * scanNx * scanNy;
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
            dfIdx, dfNum, yIdx, scanNy);
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
            WriteComplexBinaryFile(filename, wave, preferrence, varargin{:});
        end
    end
end

delete(wbHandle);

% nested functions:
    function r = CreateNewFolder(folderName)
        r = 1;
        status = mkdir(folderName);
        if status == 0
            errorMessage = sprintf('Error: creating %s failed', folderName);
            uiwait(warndlg(errorMessage));
            r = 0;
        end
    end

    function r = DirectoryExist
        r = 1;
        if ~isfolder(destDir)
            errorMessage = sprintf('Error: %s does not exist!\n', destDir);
            uiwait(warndlg(errorMessage));
            r = 0;
        end
    end

    function r = ValidAberrationType
        r = 1;
        if ~(strcmp(aberrType, 'reduced') || strcmp(aberrType, 'full'))
            errorMessage = 'Error: invalid input of aberration type!';
            uiwait(warndlg(errorMessage));
            r = 0;
        end
    end

    function r = ValidPreferrence
        r = 1;
        if ~(strcmp(preferrence, 'column') || strcmp(preferrence, 'row'))
            errorMessage = 'Error: invalid input of preferrence!';
            uiwait(warndlg(errorMessage));
            r = 0;
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

