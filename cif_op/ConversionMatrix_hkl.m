function [convMat] = ConversionMatrix_hkl(cellLengths, cellAngles, hkl)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if all(cellLengths) && all(cellAngles) && any(hkl)
    a = cellLengths(1);
    b = cellLengths(2);
    c = cellLengths(3);

    alpha = cellAngles(1);
    beta = cellAngles(2);
    gamma = cellAngles(3);

    c1 = c * cosd(alpha);
    c2 = c * (cosd(beta) - cosd(gamma) * cosd(alpha)) / sind(gamma);
    c3 = sqrt(c^2 - c1^2 - c2^2);

    initConvMat = [a, b * cosd(gamma), c1;
        0, b * sind(gamma), c2;
        0, 0, c3];

    viewDirection = hkl(1) * initConvMat(:, 1) +...
        hkl(2) * initConvMat(:, 2) +...
        hkl(3) * initConvMat(:, 3);

    rotMat = RotationOperator(viewDirection);
    convMat = rotMat * initConvMat;
else
    convMat = zeros(3);
end

end

