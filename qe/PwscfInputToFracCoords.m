function [fracCoords, a1, a2, a3, b1, b2, b3] = PwscfInputToFracCoords(pwscf)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

bohrRadius = 5.29177210903e-1;

if isfield(pwscf.system, 'celldm_1')
    % cell parameters defined with celldm in Bohr
    celldm = zeros(1, 6);
    celldm(1) = pwscf.system.celldm_1 * bohrRadius;
    alat = celldm(1);
    if isfield(pwscf.system, 'celldm_2')
        celldm(2) = pwscf.system.celldm_2;
    end
    if isfield(pwscf.system, 'celldm_3')
        celldm(3) = pwscf.system.celldm_3;
    end
    if isfield(pwscf.system, 'celldm_4')
        celldm(4) = pwscf.system.celldm_4;
    end
    if isfield(pwscf.system, 'celldm_5')
        celldm(5) = pwscf.system.celldm_5;
    end
    if isfield(pwscf.system, 'celldm_6')
        celldm(6) = pwscf.system.celldm_6;
    end

    [a1, a2, a3, b1, b2, b3] = PwscfBravisLatticeVector(pwscf.system.ibrav, celldm);
elseif isfield(pwscf.system, 'A')
    a = pwscf.system.A;
    alat = a;
    b = a;
    c = a;
    cosab = 0;
    cosac = 0;
    cosbc = 0;
    if isfield(pwscf.system, 'B')
        b = pwscf.system.B;
    end
    if isfield(pwscf.system, 'C')
        c = pwscf.system.C;
    end
    if isfield(pwscf.system, 'cosAB')
        cosab = pwscf.system.cosAB;
    end
    if isfield(pwscf.system, 'cosAC')
        cosac = pwscf.system.cosAC;
    end
    if isfield(pwscf.system, 'cosBC')
        cosbc = pwscf.system.cosBC;
    end

    [a1, a2, a3, b1, b2, b3] = PwscfBravisLatticeVector(pwscf.system.ibrav, a, b, c, cosab, cosac, cosbc);
else
    error("Invalid input");
end

a1 = reshape(a1, 1, []);
a2 = reshape(a2, 1, []);
a3 = reshape(a3, 1, []);

if strcmp(pwscf.atomic_positions.option, 'alat')
    xs = reshape(pwscf.atomic_positions.xs * alat, 1, []);
    ys = reshape(pwscf.atomic_positions.ys * alat, 1, []);
    zs = reshape(pwscf.atomic_positions.zs * alat, 1, []);
    xyz = [xs; ys; zs];
    convMat = BasesToConvMat(a1, a2, a3);
    fracXyz = convMat \ xyz;
    fracXyz = mod(fracXyz, 1);
elseif strcmp(pwscf.atomic_positions.option, 'bohr')
    xs = reshape(pwscf.atomic_positions.xs * bohrRadius, 1, []);
    ys = reshape(pwscf.atomic_positions.ys * bohrRadius, 1, []);
    zs = reshape(pwscf.atomic_positions.zs * bohrRadius, 1, []);
    xyz = [xs; ys; zs];
    convMat = BasesToConvMat(a1, a2, a3);
    fracXyz = convMat \ xyz;
    fracXyz = mod(fracXyz, 1);
elseif strcmp(pwscf.atomic_positions.option, 'angstrom')
    xs = reshape(pwscf.atomic_positions.xs, 1, []);
    ys = reshape(pwscf.atomic_positions.ys, 1, []);
    zs = reshape(pwscf.atomic_positions.zs, 1, []);
    xyz = [xs; ys; zs];
    convMat = BasesToConvMat(a1, a2, a3);
    fracXyz = convMat \ xyz;
    fracXyz = mod(fracXyz, 1);
elseif strcmp(pwscf.atomic_positions.option, 'crystal') ||...
        strcmp(pwscf.atomic_positions.option, 'crystalsg')
    fracXyz = [reshape(pwscf.atomic_positions.xs, 1, []);...
               reshape(pwscf.atomic_positions.ys, 1, []);...
               reshape(pwscf.atomic_positions.zs, 1, [])];
end

fracCoords = [reshape(pwscf.atomic_positions.types, 1, []);...
           ones(1, pwscf.system.nat);...
           fracXyz];

end