function [exitWave] = MsbtMultislice(inciWave, wavLen, fxMesh, fyMesh, a3p, ...
    transFuncs, sliceDists, nStack)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

ku = a3p / norm(a3p);

tmpWave = fftshift(inciWave);
nSlicePerStack = length(sliceDists);
[n2, n1] = size(inciWave);
shiftedPropKernels = 1i * ones(n2, n1, nSlicePerStack);
for iSlice = 1 : nSlicePerStack
    shiftedPropKernels(:, :, iSlice) = fftshift(MsbtPropKernel(fxMesh, fyMesh, ...
        ku, wavLen, sliceDists(iSlice)));
end

shiftedTransFuncs = fftshift(transFuncs, 1);
shiftedTransFuncs = fftshift(shiftedTransFuncs, 2);

for iStack = 1 : nStack
    for iSlice = 1 : nSlicePerStack
        tmpWave = tmpWave .* shiftedTransFuncs(:, :, iSlice);
        tmpWave = ifft2(shiftedPropKernels(:, :, iSlice) .* fft2(tmpWave));
    end
end

exitWave = ifftshift(tmpWave);

end