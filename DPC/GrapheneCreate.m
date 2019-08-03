function [Vz] = GrapheneCreate()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a0 = 1.42;% in Angstrom
a = a0 * sqrt(3);
% Graphene structure on XY plane
LAx = a0 * [0 1 3 1.5 2.5 0 1 3];
LAy = a * [0 0 0 0.5 0.5 1 1 1];
LA1 = [LAx; LAy];
LA2 = [LAx - 3 * a0 * ones(1, length(LAx)); LAy];
LA3 = [LAx - 3 * a0 * ones(1, length(LAx)); LAy - a * ones(1, length(LAy))];
LA4 = [LAx; LAy - a * ones(1, length(LAy))];

% Extend these hexagons to fully cover the scope
% Extend the first two slices, say A and B, along x axis
LA = [LA1 LA2 LA3 LA4];

LA0 = [LA1 LA4];  LA0p=[LA2 LA3];
for i = 1 : 2
    LA = [LA LA0 + [i * 3 * a0 * ones(1, size(LA0, 2)); zeros(1, size(LA0, 2))] LA0p + [-i * 3 * a0 * ones(1, size(LA0p, 2)); zeros(1, size(LA0p, 2))]];
end
LA0=[LA1 LA2];  LA0p=[LA3 LA4];
for i = 1 : 2
    LA0 = [LA0 LA1 + [i * 3 * a0 * ones(1, size(LA1, 2)); zeros(1, size(LA1, 2))] LA2 + [-i * 3 * a0 * ones(1, size(LA2, 2)); zeros(1, size(LA2, 2))]];
    LA0p = [LA0p LA4 + [i * 3 * a0 * ones(1, size(LA4, 2)); zeros(1, size(LA4, 2))] LA3 + [-i * 3 * a0 * ones(1, size(LA3, 2)); zeros(1, size(LA3, 2))]];
end

% Extend three slices along y axis
for i = 1 : 4
    LA = [LA LA0 + [zeros(1, size(LA0, 2)); i * a * ones(1, size(LA0, 2))] LA0p + [zeros(1, size(LA0p, 2)); -i * a * ones(1, size(LA0p, 2))]];
end
% Delete the extra elements
LA=(uniquetol(LA.', 'ByRows', true)).';

Nx = 512;
Ny = 512;
Lx = 6 * 3 * a0;
Ly = 10 * a;
Vz = zeros(Ny, Nx);
for i=1:size(LA,2)
    Vz = Vz + ProjectedPotential(Lx, Ly, Nx, Ny, 6, LA(1,i), LA(2,i));
end

end

