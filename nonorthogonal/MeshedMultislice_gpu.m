function exitWave = MeshedMultislice_gpu(inciWave, wavLen, fxMesh, fyMesh, transFuncs, ...
    sliceDists, nStack)

tmpWave = fftshift(inciWave);
nSlicePerStack = length(sliceDists);
[n2, n1] = size(inciWave);
shiftedPropKernels = 1i * ones(n2, n1, nSlicePerStack, "gpuArray");

for iSlice = 1 : nSlicePerStack
    shiftedPropKernels(:, :, iSlice) = fftshift(MeshedFresnelPropKernel(fxMesh, ...
        fyMesh, wavLen, sliceDists(iSlice)));
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