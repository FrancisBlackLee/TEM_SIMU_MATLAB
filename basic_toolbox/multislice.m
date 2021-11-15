function [varargout] = multislice(incidentWave, wavLen, Lx, Ly,...
    transFuncs, sliceDist, stackNum, saveSliceSeries, saveDir)
%MULTISLICE.M performs the multislice procedure. See E. J. Kirkland
%Advanced Computing in Electron Microscopy.
% Input:
%   IncidentWave;
%   Lx, Ly -- sampling parameters;
%   TransFuncs -- transmission functions;
%   LayerDist -- distance from one layer to the next layer;
%   StackNum -- number of stackings;
%   SaveSliceSeries -- Slice indices upon which the wave functions are to
%       be saved, i.e. [1, 2, 6, 8, 12, 34], note that the array should be
%       in an increasing order and it maximum element should be no greater
%       than N + 1, where N = StackNum * LayerNum;
%   SaveDir -- the directory where you want these wave functions to be
%       saved;
%
% Output:
%   [wave] = multislice(...) -- (1) output the exit wave if SaveSliceSeries
%       is not specified; (2) concatenate the selected wave functions to
%       the exit wave function if SaveSliceSeries is specified but saving
%       data to saveDir fails or saveDir is not specified;
%   [wave, selectWaves] -- output the exit wave and the selected wave
%       functions if SaveSliceSeries is specified.


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

bwlProp = 2.0/3.0;

switch nargin
    case 7
        wave = incidentWave;
        sliceNum = length(sliceDist);
        for stackIdx = 1: stackNum
            for sliceIdxInStack = 1: sliceNum
                wave = wave .* transFuncs(:, :, sliceIdxInStack);
                wave = propTF_1(wave, Lx, Ly, wavLen, sliceDist(sliceIdxInStack), bwlProp);
            end
        end
        wave = wave .* transFuncs(:, :, 1);
        varargout{1} = wave;
    otherwise
        if nargin == 8
            saveWaveMat = 0;
        elseif nargin == 9
            saveWaveMat = 1;
            mkDirStat = mkdir(saveDir);
            if (mkDirStat == 0) && (~isfolder(saveDir))
                saveWaveMat = 0;
                warning(['Creating output folder failed!',...
                    ' Selected slice wave functions are concatenated to the output argument']);
            end
        end
        
        saveSliceIdx = 1;
        saveSliceNum = length(saveSliceSeries);
        wave = incidentWave;
        [Ny, Nx] = size(wave);
        waveMat_3d = 1i * ones(Ny, Nx, saveSliceNum);
        sliceNum = length(sliceDist);
        for stackIdx = 1: stackNum
            for sliceIdxInStack = 1: sliceNum
                wave = wave .* transFuncs(:, :, sliceIdxInStack);
                wave = propTF_1(wave, Lx, Ly, wavLen, sliceDist(sliceIdxInStack), bwlProp);
                
                sliceIdx = (stackIdx - 1) * sliceNum + sliceIdxInStack;
                if saveSliceIdx <= saveSliceNum
                    if sliceIdx == saveSliceSeries(saveSliceIdx)
                        waveMat_3d(:, :, saveSliceIdx) = wave;
                        saveSliceIdx = saveSliceIdx + 1;
                    end
                end
            end
        end
        
        if saveWaveMat == 1
            filename = fullfile(saveDir, 'wave_data.mat');
            save(filename, 'waveMat_3d', 'saveSliceSeries', '-v7.3');
        end
        
        if nargout == 1
            if saveWaveMat == 1
                varargout{1} = wave;
            else
                varargout{1} = cat(3, wave, waveMat_3d);
            end
        elseif nargout == 2
            varargout{1} = wave;
            varargout{2} = waveMat_3d;
        end
end

end

