function [convMat, glideVec] = SymmetryOperatorToMatrix(stringExpression)
%SymmetryOperatorToMatrix() calculates the symmetry operator as a matrix by
%parsing the input 3 expressions, separated by commas.
% Input:
%   stringExpression -- symmetry operator as 3 expressions;
% Output:
%   matrixExpression -- symmetry operator as a 3-by-3 matrix;

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

