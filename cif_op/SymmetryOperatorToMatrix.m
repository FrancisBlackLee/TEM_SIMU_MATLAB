function [convMat, glideVec] = SymmetryOperatorToMatrix(stringExpression)
%SymmetryOperatorToMatrix() calculates the symmetry operator as a matrix by
%parsing the input 3 expressions, separated by commas.
% Input:
%   stringExpression -- symmetry operator as 3 expressions;
% Output:
%   matrixExpression -- symmetry operator as a 3-by-3 matrix;

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

convMat = zeros(3);
glideVec = zeros(3, 1);

% find commas in the expression and split the expression
commaIndices = strfind(stringExpression, ',');
rowStr_1 = stringExpression(1 : commaIndices(1) - 1);
rowStr_2 = stringExpression(commaIndices(1) + 1 : commaIndices(2) - 1);
rowStr_3 = stringExpression(commaIndices(2) + 1 : end);

% process rowStr_1, rowStr_2 and rowStr_3
[convMat(1, :), glideVec(1)] = SymmetryOperatorToRowVec(rowStr_1);
[convMat(2, :), glideVec(2)] = SymmetryOperatorToRowVec(rowStr_2);
[convMat(3, :), glideVec(3)] = SymmetryOperatorToRowVec(rowStr_3);

end

