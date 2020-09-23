function [rowVec, glide] = SymmetryOperatorToRowVec(rowStr)
%SymmetryOperatorToRowVec() calculates the symmetry operator as a vector by
%parsing the input expression.
% Input:
%   rowStr -- symmetry operator an expression;
% Output:
%   rowVec -- symmetry operator as a 1-by-3 vector;

rowVecX = [1, 0, 0];
rowVecY = [0, 1, 0];
rowVecZ = [0, 0, 1];

rowVec = zeros(1, 3);
% find x and its sign:
if contains(rowStr, 'x')
    rowStr_xIdx = strfind(rowStr, 'x');
    if rowStr_xIdx > 1
        if rowStr(rowStr_xIdx - 1) == '-'
            rowVec = rowVec - rowVecX;
        else
            rowVec = rowVec + rowVecX;
        end
        rowStr(rowStr_xIdx - 1 : rowStr_xIdx) = [];
    else
        rowVec = rowVec + rowVecX;
        rowStr(rowStr_xIdx) = [];
    end
end

% find y and its sign:
if contains(rowStr, 'y')
    rowStr_yIdx = strfind(rowStr, 'y');
    if rowStr_yIdx > 1
        if rowStr(rowStr_yIdx - 1) == '-'
            rowVec = rowVec - rowVecY;
        else
            rowVec = rowVec + rowVecY;
        end
        rowStr(rowStr_yIdx - 1 : rowStr_yIdx) = [];
    else
        rowVec = rowVec + rowVecY;
        rowStr(rowStr_yIdx) = [];
    end
end

% find z and its sign:
if contains(rowStr, 'z')
    rowStr_zIdx = strfind(rowStr, 'z');
    if rowStr_zIdx > 1
        if rowStr(rowStr_zIdx - 1) == '-'
            rowVec = rowVec - rowVecZ;
        else
            rowVec = rowVec + rowVecZ;
        end
        rowStr(rowStr_zIdx - 1 : rowStr_zIdx) = [];
    else
        rowVec = rowVec + rowVecZ;
        rowStr(rowStr_zIdx) = [];
    end
end

if ~isempty(rowStr)
    % only glide is left, find its sign:
    glideSign = 1;
    if contains(rowStr, '+')
        plusIdx = strfind(rowStr, '+');
        rowStr(plusIdx) = [];
    elseif contains(rowStr, '-')
        glideSign = -1;
        minusIdx = strfind(rowStr, '-');
        rowStr(minusIdx) = [];
    end

    fracCell = regexp(rowStr, '[0-9]+', 'match');
    fracNumerator = str2double(fracCell{1});
    fracDenominator = str2double(fracCell{2});

    glide = fracNumerator / fracDenominator;
else
    glide = 0;
end

end

