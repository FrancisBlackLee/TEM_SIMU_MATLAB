function [sliceSeImgs] = STEM_X_slice_SE_var_z_gpu(Lx, Ly, params, transFuncs, seObjFuncs,...
    mfp, sliceDist, stackNum, z0, aberrType)
% output se contribution of each slice

if nargin == 9
    aberrType = 'reduced';
end
[Ny, Nx, sliceNumPerStack] = size(transFuncs);
dx = Lx / Nx;
dy = Ly / Ny;
wavLen = HighEnergyWavLen_X(params.KeV);
% generate fftshifted Fresnel propagation kernels:
shiftPropKer = 1i * ones(Ny, Nx, sliceNumPerStack, "single", "gpuArray");
for sliceIdx = 1 : sliceNumPerStack
    shiftPropKer(:, :, sliceIdx) = fftshift(gpuArray(single(FresnelPropKernel_X(Lx, Ly,...
        Nx, Ny, wavLen, sliceDist(sliceIdx)))));
end
% generate objective transfer function with an aperture:
if strcmp(aberrType, 'reduced')
    otf = params.aperture .* ObjTransFunc_X(params, Lx, Ly, Nx, Ny);
elseif strcmp(aberrType, 'full')
    otfPhase = AberrationPhaseShift_X(params.aberration, wavLen, Lx, Ly, Nx, Ny);
    otf = params.aperture .* exp(-1i * otfPhase);
else
    error('Invalid aberration type!\n');
end

% fftshift all the transmission function in place:
transFuncs = fftshift(transFuncs, 1);
transFuncs = fftshift(transFuncs, 2);

% fftshift all the secondary electron object functions in place:
seObjFuncs = fftshift(seObjFuncs, 1);
seObjFuncs = fftshift(seObjFuncs, 2);

% start scanning:
scanNx = length(params.scanx);
scanNy = length(params.scany);

totalSliceNum = stackNum * sliceNumPerStack;
sliceSeImgs = zeros(scanNy, scanNx, totalSliceNum, "single", "gpuArray");

process = waitbar(0, 'start scanning');
for iy = 1 : scanNy
    for ix = 1 : scanNx
        tempWave = fftshift(GenerateProbe_X(otf, params.scanx(ix), params.scany(iy),...
            Lx, Ly, Nx, Ny));

        depth = 0.0;
        for stackIdx = 1 : stackNum
            for sliceIdx = 1 : sliceNumPerStack
                totalSliceIdx = (stackIdx - 1) * sliceNumPerStack + sliceIdx;
                tempWave = tempWave .* transFuncs(:,:,sliceIdx);
                tempWave = ifft2(shiftPropKer(:, :, sliceIdx) .* fft2(tempWave));
                depth = depth + sliceDist(sliceIdx);
                
                generationRate = sum((abs(tempWave.^2) .* seObjFuncs(:, :, sliceIdx)), 'all');
                seTransport = exp(-(depth - z0(iy, ix)) / mfp(iy, ix));
                seSignal = generationRate * seTransport;
                sliceSeImgs(iy, ix, totalSliceIdx) = seSignal;
            end
        end
    end

    doneRatio = iy / scanNy;
    waitbar(doneRatio, process, [num2str(roundn(doneRatio, -3) * 100), '%']);
end
delete(process);

end