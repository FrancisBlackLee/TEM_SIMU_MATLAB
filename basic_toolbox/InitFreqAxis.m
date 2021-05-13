function [f] = InitFreqAxis(L, N)
%InitFreqAxis.m initializes the spatial frequency axis with the origin at
%its center.
% Input:
%   L -- side length;
%   N -- sampling number;
% Output:
%   f -- spatial frequency axis;

d = L / N;
f = -1 / (2 * d) : 1 / L : 1 / (2 * d) - 1 / L;

end

