function [varargout] = ConvMatToBases(convMat)
%ConvMatToBases.m converts the conversion matrix to the lattice or cell
%bases.
% Input:
%   convMat -- conversion matrix;
% Output:
%   bases.a -- reoriented lattice vector A in cartesian coordinates;
%   bases.b -- reoriented lattice vector B in cartesian coordinates;
%   bases.c -- reoriented lattice vector C in cartesian coordinates;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2023  Francis Black Lee (Li Xian)

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

if nargout == 1
    bases.a = convMat(:, 1)';
    bases.b = convMat(:, 2)';
    bases.c = convMat(:, 3)';
    varargout{1} = bases;
elseif nargout == 3
    a1 = convMat(:, 1)';
    a2 = convMat(:, 2)';
    a3 = convMat(:, 3)';
    varargout{1} = a1;
    varargout{2} = a2;
    varargout{3} = a3;
else
    error('Incorrect number of output arguments');
end

end


