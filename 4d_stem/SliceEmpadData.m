function SliceEmpadData(srcFilename, destDir, tNum, zNum,...
    destFilenamePrefix, destFilenameIdxDigit)
%SliceEmpadData.m slices the EMPAD data into single images (CBED patterns)
%with respect to each scanning position.
%   srcFilename -- filename of EMPAD data;
%   destDir -- destination directory where the sliced images will be saved.
%   tNum, zNum -- number of scanning positions per dimension, t: dimension
%       1 (vertical), z: dimension 2 (horizontal);
%   destFilenamePrefix -- filename prefix of the sliced images;
%   destFilenameIdxDigit -- number of digits for indexing the scanning
%       positions in the filename

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

[fileID, errMsg] = fopen(srcFilename);
if fileID == -1
    error(errMsg);
end

imageSize = [128, 130];
cropRangeY = 3 : 126;
cropRangeX = 3 : 126;

idxFormat = ['%0', num2str(destFilenameIdxDigit), 'd'];
for tIdx = 1 : tNum
    for zIdx = 1 : zNum
        tmpRawImage = fread(fileID, imageSize, 'float', 'ieee-le');
        tmpCropImage = tmpRawImage(cropRangeY, cropRangeX)';
        destFilename = [destFilenamePrefix, '_t', num2str(tIdx, idxFormat),...
            '_z', num2str(zIdx, idxFormat), '.raw'];
        destFilename = fullfile(destDir, destFilename);
        WriteBinaryFile(destFilename, tmpCropImage, 'row', 'float', 'ieee-be');
    end
end

fclose(fileID);

end

