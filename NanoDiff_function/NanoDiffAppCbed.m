function [diffPatt] = NanoDiffAppCbed(keV, objTransFunc, probeX, probeY,...
    projPotMat, sliceDist, pixelSize, dimPixelNum)
%NanoDiffAppCbed() computes the convergent beam electron diffraction
%pattern as an external API for NanoDiff app.
% Input:
%   keV -- electron beam energy;
%   objTransFunc -- objective transfer function (already including
%       aperture);
%   xp, yp
%   projPotMat -- projected potential matrix;
%   sliceDist -- slice distance list;
% Output:
%   diffPatt -- diffraction pattern;

sideLength = pixelSize * dimPixelNum;
interCoeff = InteractionCoefficient(keV);
sliceNum = length(sliceDist);

transFuncMat = 1i * ones(dimPixelNum, dimPixelNum, sliceNum);
for sliceIdx = 1 : sliceNum
    transFuncMat(:, :, sliceIdx) =...
        exp(1i * interCoeff * projPotMat(:, :, sliceIdx) / 1e3);
end

wave = GenerateProbe_X(objTransFunc, probeX, probeY, sideLength, sideLength,...
    dimPixelNum, dimPixelNum);

stackNum = 1;
wave = multislice_X(wave, keV, sideLength, sideLength, transFuncMat,...
    sliceDist, stackNum);

wave = ifftshift(fft2(fftshift(wave)));
diffPatt = abs(wave.^2);
diffPatt = log(1 + 0.1 * diffPatt);

end

