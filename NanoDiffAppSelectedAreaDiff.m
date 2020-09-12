function [diffPatt] = NanoDiffAppSelectedAreaDiff(keV, projPotMat,...
    sliceDist, pixelSize, dimPixelNum)
%NanoDiffAppSelectedAreaDiff() computes the selected area diffraction
%pattern as an external API for NanoDiff app.
% Input:
%   keV -- electron beam energy;
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

wave = ones(dimPixelNum, dimPixelNum);
stackNum = 1;
wave = multislice_X(wave, keV, sideLength, sideLength, transFuncMat,...
    sliceDist, stackNum);

wave = ifftshift(fft2(fftshift(wave))) * pixelSize^2;
diffPatt = abs(wave.^2);
diffPatt = log(1 + 0.1 * diffPatt);

end

