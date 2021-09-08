function [data] = ReadComplexBinaryFile(filename, dataSize, preferrence, varargin)
%ReadComplexBinaryFile.m read a complex matrix from a binary file, in which
%the real components start from the first element with stride = 2, and the
%imaginary components start from the second element with stride = 2, in
%accord with the protocol of vtemlab. Note: only 2-dimension matrices are
%valid, the output data is converted to matlab complex matrix.
%   filename -- Binary file name;
%   dataSize -- data size;
%   preferrence -- writing by column (default) or row;
%   varargin -- specific parameters to define the reading mode (default:
%       double);

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

Ny = dataSize(1);
Nx = dataSize(2);
tmpDataSize = [Ny, 2 * Nx];

fileID = fopen(filename);
if nargin == 2
    tmpData = fread(fileID, tmpDataSize, 'double');
elseif nargin > 2
    if strcmp(preferrence, 'column')
        % do nothing
    elseif strcmp(preferrence, 'row')
        tmpDataSize = [2 * Nx, Ny];
    else
        disp('Error: invalid input of preferrence!');
        return;
    end
    
    if nargin == 3
        tmpData = fread(fileID, tmpDataSize, 'double');
    else
        tmpData = fread(fileID, tmpDataSize, varargin{:});
    end
    
    if strcmp(preferrence, 'row')
        tmpData = tmpData';
    end
end
fclose(fileID);

data = tmpData(:, 1 : 2 : 2 * Nx - 1) + 1i * tmpData(:, 2 : 2 : 2 * Nx);

end

