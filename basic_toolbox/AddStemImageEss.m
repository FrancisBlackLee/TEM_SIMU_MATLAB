function [essImg] = AddStemImageEss(rawImg, scanLx, scanLy, scanNx, scanNy, ess)
%AddStemImageEss adds effective source size to stem image.
%   To be added.

scanDx = scanLx / (scanNx - 1);
scanDy = scanLy / (scanNy - 1);
sigma = ess / 2.35482 / sqrt(scanDx * scanDy);
essImg = imgaussfilt(rawImg, sigma, 'Padding', 'symmetric');

end