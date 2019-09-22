function [wasted_ef] = wasted_temp()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
L = 20;     %side length 
M = 256;      %number of samples
dx = L/(M-1);   %scr sample interval
x = linspace(-L/2,L/2,M);  %scr coords
y = x;
Proj_Pot = ProjectedPotential(L, L, M, M, 6, -8, 0) + ...
           ProjectedPotential(L, L, M, M, 14, -4, 0) + ...
           ProjectedPotential(L, L, M, M, 29, 0, 0) + ...
           ProjectedPotential(L, L, M, M, 79, 4, 0) + ...
           ProjectedPotential(L, L, M, M, 92, 8, 0);
[EleField_x, EleField_y] = gradient(Proj_Pot/1000, dx, dx);
wasted_ef = sqrt(EleField_x.^2 + EleField_y.^2);

end

