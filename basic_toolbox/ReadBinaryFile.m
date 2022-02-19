function [data] = ReadBinaryFile(filename, dataSize, preferrence, varargin)
%ReadBinaryFile reads matrix from the destination binary file.
%   filename -- Binary file name;
%   dataSize -- data size;
%   preferrence -- writing by column (default) or row;
%   varargin -- specific parameters to define the reading mode (default:
%       double);

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

fileID = fopen(filename);
if nargin == 2
    data = fread(fileID, dataSize, 'double');
elseif nargin > 2
    Ny = dataSize(1);
    Nx = dataSize(2);
    if strcmp(preferrence, 'column')
        % do nothing
    elseif strcmp(preferrence, 'row')
        dataSize = [Nx, Ny];
    else
        disp('Error: invalid input of preferrence!');
        return;
    end
    
    if nargin == 3
        data = fread(fileID, dataSize, 'double');
    else
        data = fread(fileID, dataSize, varargin{:});
    end
    
    if strcmp(preferrence, 'row')
        data = data';
    end
end
fclose(fileID);

end


