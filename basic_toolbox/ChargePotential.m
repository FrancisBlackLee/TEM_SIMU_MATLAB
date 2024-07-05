function [p] = ChargePotential(q, xq, yq, zq, x, y, z, thr)
%ChargePotential() calculates electric potential of charges
%   q -- charges (C)
%   xq, yp, zq -- position of the charges (angstrom)
%   x, y, z -- coordinate grids (angstrom)
%   thr -- threshold for removing the singularity

p = zeros(size(x));
nq = length(q);
for iq = 1 : nq
    r = sqrt((x - xq(iq)).^2 + (y - yq(iq)).^2 + (z - zq(iq)).^2);
    r(r < thr) = thr;
    p = p + q(iq) ./ r;
end

epsilon0 = 8.8541878188e-22; % F / angstrom
p = p  / (4 * pi * epsilon0);

end