function [convMat] = BasesToConvMat(varargin)
%BasesToConvMat converts the lattice or cell bases to the conversion
%matrix.
% [convMat] = BasesToConvMat(bases):
%   bases.a -- lattice vector A in cartesian coordinates;
%   bases.b -- lattice vector B in cartesian coordinates;
%   bases.c -- lattice vector C in cartesian coordinates;
%
% [convMat] = BasesToConvMat(a1, a2, a3):
%   a1 -- lattice vector A in cartesian coordinates;
%   a2 -- lattice vector B in cartesian coordinates;
%   a3 -- lattice vector C in cartesian coordinates;
%

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

if nargin == 1
    bases = varargin{1};
    vecA = reshape(bases.a, [], 1);
    vecB = reshape(bases.b, [], 1);
    vecC = reshape(bases.c, [], 1);
elseif nargin == 3
    vecA = reshape(varargin{1}, [], 1);
    vecB = reshape(varargin{2}, [], 1);
    vecC = reshape(varargin{3}, [], 1);
else
    error('Incorrect number of input arguments');
end

convMat = [vecA, vecB, vecC];

end


