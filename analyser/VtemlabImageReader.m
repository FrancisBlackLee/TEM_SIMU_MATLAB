function [imageMatrix] = VtemlabImageReader(filename, size1, size2, varargin)
%VtemlabImageReader.m converts the vtemlab image data to matlab matrix.
%   varargin -- extra arguments for reading binary file except the
%       precision option.

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

[~, ~, ext] = fileparts(filename);
if strcmp(ext, '.txt')
    rawData = load(filename);
    size3 = size(rawData, 1);
    imageMatrix = zeros(size1, size2, size3);
    for imageIdx = 1 : size3
        tmpImage = reshape(rawData(imageIdx, :), size2, size1);
        imageMatrix(:, :, imageIdx) = tmpImage';
    end
    
elseif strcmp(ext, '.bin')
    fileID = fopen(filename);
    imageNum = fread(fileID, 1, 'int', varargin{:});
    rowNum = fread(fileID, 1, 'int', varargin{:});
    colNum = fread(fileID, 1, 'int', varargin{:});
    imageMatrix = zeros(rowNum, colNum, imageNum);
    for imageIdx = 1 : imageNum
        tmpImage = fread(fileID, [colNum, rowNum], 'double', varargin{:});
        imageMatrix = tmpImage';
    end
end

end

