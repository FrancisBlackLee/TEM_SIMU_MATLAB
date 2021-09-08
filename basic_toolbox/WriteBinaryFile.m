function WriteBinaryFile(filename, data, preferrence, varargin)
%WriteBinaryFile writes a matrix as binary files to the destination place.
%   filename -- Binary file name;
%   data -- matrix to be saved.
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

fileID = fopen(filename, 'w');
if nargin == 2
    fwrite(fileID, data, 'double');
elseif nargin > 2
    if strcmp(preferrence, 'column')
        % do nothing
    elseif strcmp(preferrence, 'row')
        data = data';
    else
        disp('Error: invalid input of preferrence!');
        return;
    end
    
    if nargin == 3
        fwrite(fileID, data, 'double');
    else
        fwrite(fileID, data, varargin{:});
    end
end
fclose(fileID);

end

