function [x] = InitAxis(L, N)
%InitAxis.m initializes the spatial axis with the origin at its center.
% Input:
%   L -- side length;
%   N -- sampling number;
% Output:
%   x -- spatial axis;

d = L / N;
x = -L / 2 : d : L / 2 - d;

end

