function [exitWave] = multislice_X(inciWave, KeV, Lx, Ly, transFuncs,...
    sliceDist, stackNum, projPotDir, fileExtension)
%multislice_X.m performs the multislice procedure. See E. J. Kirkland
%Advanced Computing in Electron Microscopy for more details.
%   InciWave -- incident wave;
%   KeV -- Energy of incident electron beam (in KeV);
%   Lx, Ly -- sampling side lengths;
%   TransFuncs -- 3D array including transmission functions,
%       TransFuncs(:, :, i) denotes the ith transmission function. Note
%       that for a bulk (large volume) material possibly containing too
%       many transmission functions, creating such a 3D array might cause
%       memory overflow, thus input TransFuncs = 'files', the program will
%       load projected potential files under the given path (the optional
%       input), considering that for such bulk materials, these slices do
%       not need to be looped, these files are only loaded once.
%   SliceDist -- array whose elements are distances from the identically
%       indexed slice to the next slice;
%   StackNum -- number of stackings for the slices;
%   ProjPotDir -- directory where the projected potential (in V-Angs) are
%       stored. Also note that these ProjPot files are named in the same
%       style, name ordering is the same as the slice ordering.
%   FileExtension -- a required input if TransFuncDir is input. '*.txt' is
%       suggested.
%   Note: X denotes an experimental version. This function can be further
%       optimized for ADF-STEM, but there is no need to keep modifying this
%       function, thus a dependent multislice for ADF-STEM is to be
%       released.

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

wavLen = HighEnergyWavLen_X(KeV);
switch nargin
    case 7
        tempWave = fftshift(inciWave);
        sliceNum = length(sliceDist);
        [Ny, Nx] = size(inciWave);
        shiftedPropKernels = zeros(Ny, Nx, sliceNum) + 1i * ones(Ny, Nx, sliceNum);
        for sliceIdx = 1 : sliceNum
            shiftedPropKernels(:, :, sliceIdx) = ...
                fftshift(FresnelPropKernel_X(Lx, Ly, Nx, Ny, wavLen, sliceDist(sliceIdx)));
        end
        for stackIdx = 1 : stackNum
            for sliceIdx = 1 : sliceNum
                tempWave = tempWave .* fftshift(transFuncs(:, :, sliceIdx));
                tempWave = ifft2(shiftedPropKernels(:, :, sliceIdx) .* fft2(tempWave));
            end
        end
        exitWave = ifftshift(tempWave);
    otherwise
        if ~isfolder(projPotDir)
          errorMessage = sprintf('Error: The following folder does not exist:\n%s', projPotDir);
          uiwait(warndlg(errorMessage));
          return;
        end
        interCoeff = InteractionCoefficient(KeV);
        projPotFiles = dir(fullfile(projPotDir,fileExtension));
        filenames = {projPotFiles.name}';
        sortedNames = natsortfiles(filenames);
        tempWave = fftshift(inciWave);
        [Ny, Nx] = size(inciWave);
        for fileIdx = 1 : min(numel(sortedNames), length(sliceDist))
            filename = fullfile(projPotDir, sortedNames{fileIdx});
            if strcmp(fileExtension, '*.txt')
                tempProjPot = load(filename);
            elseif strcmp(fileExtension, '*.bin')
                fileID = fopen(filename);
                tempProjPot = fread(fileID, [Ny, Nx], 'double');
                fclose(fileID);
            else
                tempProjPot = zeros(Ny, Nx);
            end
            
            tempTransFunc = exp(1i * interCoeff * tempProjPot / 1e3);
            tempWave = tempWave .* fftshift(tempTransFunc);
            shiftedPropKernel = fftshift(FresnelPropKernel_X(Lx, Ly, Nx, Ny,...
                wavLen, sliceDist(fileIdx)));
            tempWave = ifft2(shiftedPropKernel .* fft2(tempWave));
            disp(filename);
        end
        exitWave = ifftshift(tempWave);
end

end

