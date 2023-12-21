function [superCell] = TileAsymUnitCell(unitCell, tiles, useVtemlabStyle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 2
    useVtemlabStyle = false;
end

bases = ConvMatToBases(eye(3));

na = tiles(1);
nb = tiles(2);
nc = tiles(3);

cellNum = na * nb * nc;
atomNumPerUnitCell = size(unitCell, 2);
superCell = zeros(5, cellNum * atomNumPerUnitCell);

if any(mod(tiles, 1))
    error('Invalid tiles');
end

% tile
if useVtemlabStyle
    % tile along a
    for ta = 0 : na - 1
        tmpCell = unitCell;
        tmpCell(3, :) = tmpCell(3, :) + ta;
        head = ta * atomNumPerUnitCell + 1;
        rear = head + atomNumPerUnitCell - 1;
        superCell(:, head : rear) = tmpCell;
    end

    % tile along b
    baseCell = superCell(:, 1 : na * atomNumPerUnitCell);
    for tb = 1 : nb - 1
        tmpCell = baseCell;
        tmpCell(4, :) = tmpCell(4, :) + tb;
        head = tb * na * atomNumPerUnitCell + 1;
        rear = (tb + 1) * na * atomNumPerUnitCell;
        superCell(:, head : rear) = tmpCell;
    end

    % tile along c
    baseCell = superCell(:, 1 : na * nb * atomNumPerUnitCell);
    for tc = 1 : nc - 1
        tmpCell = baseCell;
        tmpCell(5, :) = tmpCell(5, :) + tc;
        head = tc * nb * na * atomNumPerUnitCell + 1;
        rear = (tc + 1) * nb * na * atomNumPerUnitCell;
        superCell(:, head : rear) = tmpCell;
    end
else
    for tc = 0 : nc - 1
        for tb = 0 : nb - 1
            for ta = 0 : na - 1
                tmpCell = unitCell;
                tmpCell(3 : 5, :) = tmpCell(3 : 5, :) + ta * bases.a' +...
                    tb * bases.b' + tc * bases.c';
                head = (tc * na * nb + tb * na + ta) * atomNumPerUnitCell + 1;
                rear = head + atomNumPerUnitCell - 1;
                superCell(:, head : rear) = tmpCell;
            end
        end
    end
end

end