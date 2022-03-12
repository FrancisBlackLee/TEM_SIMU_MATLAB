function [wave] = DataIntensiveMultislice_X(inciWave, keV, Lx, Ly, srcDir,...
    fileOption, ppjName, sliceNum, sliceDists, stackNum, saveDepths,...
    saveOption, destDir)
%DataIntensiveMultislice_X.m is a substitute function for multislice when
%data are so large that memory overflow is possible, constant I/O is thus
%the alternative of storing all data in the workspace. Note: 1. All files
%should be  binary files; 2. Equally slicing the specimen can speed up the 
%computation.
% Input:
%   inciWave -- incident wave;
%   keV -- high tension or voltage in keV;
%   Lx, Ly -- real-space lengths for simulation;
%   srcDir -- directory of the projected potential files;
%   fileOption -- preference for reading and saving binary files:
%       fileOption.preference -- 'column'(default) / 'row';
%       fileOption.others -- extra parameters for reading or saving binary
%           files, for more details please refer to fwrite and fread;
%   ppjName -- filename prefix of the projected potential files;
%   sliceNum -- number of slices;
%   sliceDists -- distances between each slice, namely slice spacing;
%   stackNum (optional) -- number of stacks (default: 1);
%   saveDepths (optional) -- slice depths where the wave or the wave
%       intensity is saved (default: the specimen thickness);
%   saveOption (optional) -- form of data to be saved, 'wave' (default) /
%       'intensity';
%   destDir (optional) -- destination directory to save the data (default:
%       srcDir);
% Output:
%   wave -- exit wave;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2022  Francis Black Lee (Li Xian)

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

if ~isfield(fileOption, 'preference')
    fileOption.preference = 'column';
end

if ~isfield(fileOption, 'others')
    fileOption.others = {};
end

if nargin == 9
    stackNum = 1;
    saveDepths = sum(sliceDists);
    saveOption = 'wave';
    destDir = srcDir;
elseif nargin == 10
    saveDepths = stackNum * sum(sliceDists);
    saveOption = 'wave';
    destDir = srcDir;
elseif nargin == 11
    saveOption = 'wave';
    destDir = srcDir;
elseif nargin == 12
    destDir = srcDir;
end

wavLen = HighEnergyWavLen_X(keV);
interCoeff = InteractionCoefficient(keV);
[Ny, Nx] = size(inciWave);
saveSliceIndices = DepthsToSliceIndices(sliceDists, stackNum, saveDepths);

if ~any(sliceDists - sliceDists(1))
    equalSliceDist = true;
    commonPropKernel = fftshift(FresnelPropKernel_X(Lx, Ly, Nx, Ny, wavLen,...
        sliceDists(1)));
else
    equalSliceDist = false;
end

wave = fftshift(inciWave);
shiftMode.in = false;
shiftMode.out = true;
sliceCount = 0;
depthIdx = 1;
depthNum = length(saveDepths);
for stackIdx = 1 : stackNum
    for sliceIdx = 1 : sliceNum
        filename = [ppjName, num2str(sliceIdx), '.bin'];
        filename = fullfile(srcDir, filename);
        projPot = ReadBinaryFile(filename, [Ny, Nx], fileOption.preference,...
            fileOption.others{:});
        sliceTF = exp(1i * interCoeff * projPot / 1.0e3);
        sliceTF = BandwidthLimit(sliceTF, Lx, Ly, Nx, Ny, 0.67, shiftMode);
        wave = wave .* sliceTF;
        if equalSliceDist
            wave = ifft2(commonPropKernel .* fft2(wave));
        else
            slicePropKernel = fftshift(FresnelPropKernel_X(Lx, Ly, Nx, Ny,...
                wavLen, sliceDists(sliceIdx)));
            wave = ifft2(slicePropKernel .* fft2(wave));
        end
        
        sliceCount = sliceCount + 1;
        
        if depthIdx <= depthNum
            if sliceCount == saveSliceIndices(depthIdx)
                filename = [saveOption, num2str(depthIdx), '.bin'];
                filename = fullfile(destDir, filename);
                if strcmp(saveOption, 'wave')
                    WriteComplexBinaryFile(filename, ifftshift(wave),...
                        fileOption.preference, fileOption.others{:});
                elseif strcmp(saveOption, 'intensity')
                    WriteBinaryFile(filename, ifftshift(abs(wave.^2)),...
                        fileOption.preference, fileOption.others{:});
                else
                    error('Invalid input of saveOption.');
                end
                depthIdx = depthIdx + 1;
            end
        end
        
    end
end

wave = ifftshift(wave);

end

