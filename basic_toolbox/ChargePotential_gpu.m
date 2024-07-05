function [p] = ChargePotential_gpu(q, xq, yq, zq, x, y, z, thr)
%ChargePotential() calculates electric potential of charges
%   q -- charges (C)
%   xq, yp, zq -- position of the charges (angstrom)
%   x, y, z -- coordinate grids (angstrom)
%   thr -- threshold for removing the singularity

p = zeros(size(x));
nq = length(q);
[ny, nx, nz] = size(x);
for iz = 1 : nz
    x_gpu = gpuArray(x(:, :, iz));
    y_gpu = gpuArray(y(:, :, iz));
    z_gpu = gpuArray(z(:, :, iz));
    p_gpu = zeros(ny, nx, 'gpuArray');
    for iq = 1 : nq
        r_gpu = sqrt((x_gpu - xq(iq)).^2 + (y_gpu - yq(iq)).^2 + (z_gpu - zq(iq)).^2);
        r_gpu(r_gpu < thr) = thr;
        p_gpu = p_gpu + q(iq) ./ r_gpu;
    end
    p(:, :, iz) = gather(p_gpu);
end

epsilon0 = 8.8541878188e-22; % F / angstrom
p = p  / (4 * pi * epsilon0);

end