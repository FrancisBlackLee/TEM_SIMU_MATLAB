function WriteComplexBinaryFile(filename, data, preferrence, varargin)
%WriteComplexBinaryFile.m write a complex matrix as a binary file, in which
%the real components start from the first element with stride = 2, and the
%imaginary components start from the second element with stride = 2, in
%accord with the protocol of vtemlab. Note: only 2-dimension matrices are
%valid input.
%   filename -- Binary file name;
%   data -- complex matrix to be saved.
%   preferrence -- writing by column (default) or row;
%   varargin -- specific parameters to define the writing mode (default:
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

[Ny, Nx] = size(data);
tmpData = zeros(Ny, 2 * Nx);
tmpData(:, 1 : 2 : 2 * Nx - 1) = real(data);
tmpData(:, 2 : 2 : 2 * Nx) = imag(data);

fileID = fopen(filename, 'w');
if nargin == 2
    fwrite(fileID, tmpData, 'double');
elseif nargin > 2
    if strcmp(preferrence, 'column')
        % do nothing
    elseif strcmp(preferrence, 'row')
        tmpData = tmpData';
    else
        disp('Error: invalid input of preferrence!');
        return;
    end
    
    if nargin == 3
        fwrite(fileID, tmpData, 'double');
    else
        fwrite(fileID, tmpData, varargin{:});
    end
end
fclose(fileID);

end

