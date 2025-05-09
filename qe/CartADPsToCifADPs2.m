function cifAdps = CartADPsToCifADPs2(cartAdps, a1, a2, a3)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[b1, b2, b3] = DirectToReciprocal(a1, a2, a3);
convMat = BasesToConvMat(a1, a2, a3);
qMat = diag([norm(b1), norm(b2), norm(b3)] / (2 * pi));

fracAdps = zeros(size(cartAdps));

invConvMat = inv(convMat);

nAtom = size(fracAdps, 3);
for iAtom = 1 : nAtom
    fracAdps(:, :, iAtom) = invConvMat * cartAdps(:, :, iAtom) * invConvMat';
end

invQMat = inv(qMat);
cifAdps = zeros(size(fracAdps));
for iAtom = 1 : nAtom
    cifAdps(:, :, iAtom) = invQMat * fracAdps(:, :, iAtom) * invQMat';
end

end