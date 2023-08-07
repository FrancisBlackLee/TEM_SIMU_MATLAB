function [unitCell] = GlideUnitCell(unitCell, glide)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

unitCell(3, :) = unitCell(3, :) + glide(1);
unitCell(3, :) = mod(unitCell(3, :), 1);

unitCell(4, :) = unitCell(4, :) + glide(2);
unitCell(4, :) = mod(unitCell(4, :), 2);

unitCell(5, :) = unitCell(5, :) + glide(3);
unitCell(5, :) = mod(unitCell(5, :), 3);

end