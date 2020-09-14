%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2020  Francis Black Lee and Li Xian

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
function [ExitWave] = multislice(IncidentWave, WaveLength, Lx, Ly, TransFuncs, LayerDist, StackNum, SaveDir, SaveSliceSeries)
%MULTISLICE.M performs the multislice procedure. See E. J. Kirkland
%Advanced Computing in Electron Microscopy.
%   IncidentWave;
%   Lx, Ly -- sampling parameters;
%   TransFuncs -- transmission functions;
%   LayerDist -- distance from one layer to the next layer;
%   StackNum -- number of stackings;
%   SaveDir -- the directory where you want these wave functions to be
%       saved;
%   SaveSliceSeries -- Slice indices upon which the wave functions are to
%       be saved, i.e. [1, 2, 6, 8, 12, 34], note that the array should be
%       in an increasing order and it maximum element should be no greater
%       than N + 1, where N = StackNum * LayerNum;

switch nargin
    case 7
        Wave = IncidentWave;
        LayerNum = length(LayerDist);
        for i = 1: StackNum
            for j = 1: LayerNum
                Wave = Wave .* TransFuncs(:, :, j);
                Wave = propTF_1(Wave, Lx, Ly, WaveLength, LayerDist(j));
            end
        end
        ExitWave = Wave;
    otherwise
        MkDirStat = mkdir(SaveDir);
        SaveSliceIdx = 1;
        SaveSliceNum = length(SaveSliceSeries);
        Wave = IncidentWave;
        LayerNum = length(LayerDist);
        for i = 1: StackNum
            for j = 1: LayerNum
                SliceIdx = (i - 1) * LayerNum + j;
                if SaveSliceIdx <= SaveSliceNum
                    if SliceIdx == SaveSliceSeries(SaveSliceIdx)
                        FileNameReal = strcat(SaveDir, '\', num2str(SliceIdx), '_real.txt');
                        FileNameImag = strcat(SaveDir, '\', num2str(SliceIdx), '_imag.txt');
                        WaveReal = real(Wave);
                        WaveImag = imag(Wave);
                        save(FileNameReal, 'WaveReal', '-ascii', '-double', '-tabs');
                        save(FileNameImag, 'WaveImag', '-ascii', '-double', '-tabs');
                        SaveSliceIdx = SaveSliceIdx + 1;
                    end
                end
                Wave = Wave .* TransFuncs(:, :, j);
                Wave = propTF_1(Wave, Lx, Ly, WaveLength, LayerDist(j));
            end
        end
        if SaveSliceIdx <= SaveSliceNum
            if SaveSliceSeries(SaveSliceIdx) == SliceIdx + 1
                FileNameReal = strcat(SaveDir, '\', num2str(SliceIdx+1), '_real.txt');
                FileNameImag = strcat(SaveDir, '\', num2str(SliceIdx+1), '_imag.txt');
                WaveReal = real(Wave);
                WaveImag = imag(Wave);
                disp([SliceIdx + 1, SaveSliceSeries(SaveSliceIdx)]);
                save(FileNameReal, 'WaveReal', '-ascii', '-double', '-tabs');
                save(FileNameImag, 'WaveImag', '-ascii', '-double', '-tabs');
            end
        end
        ExitWave = Wave;
end

end

