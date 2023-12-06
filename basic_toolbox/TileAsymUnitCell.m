function [superCell] = TileAsymUnitCell(unitCell, tiles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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