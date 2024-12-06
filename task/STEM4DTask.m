classdef STEM4DTask
    %STEM4DTask 4D-STEM task (no TDS)
    %   Detailed explanation goes here

    properties
        keV{mustBeNumeric}
        lx{mustBeNumeric}
        ly{mustBeNumeric}
        reducedAberrations
        fullAberrations
        aperture
        scanx
        scany
        nBin{mustBeInteger}
        depths
        transFuncs
        sliceDists
        nStack{mustBeInteger}
        useFullAberration{mustBeNumericOrLogical}
        cbedDir
        useGPU{mustBeNumericOrLogical}
    end

    methods
        function obj = STEM4DTask(outDir)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.keV = 300;
            obj.reducedAberrations = InitReducedAberrations_X();
            obj.fullAberrations = InitObjectiveLensAberrations_X();
            obj.nBin = 2;
            obj.nStack = 1;
            obj.useFullAberration = true;
            obj.cbedDir = outDir;
            obj.useGPU = false;
        end

        function ExecuteTask(obj)

            if ~isfolder(obj.cbedDir)
                obj.cbedDir = 'tmpcbeddir';
                err = mkdir(obj.cbedDir);
            end

            [ny, nx, nSlicePerStack] = size(obj.transFuncs);
            wavLen = HighEnergyWavLen_X(obj.keV);
            if obj.useGPU
                shiftPropKer = 1i * ones(ny, nx, nSlicePerStack, "single", "gpuArray");
                for iSlice = 1 : nSlicePerStack
                    shiftPropKer(:, :, iSlice) = fftshift(gpuArray(single(FresnelPropKernel_X(obj.lx, ...
                        obj.ly, nx, ny, wavLen, obj.sliceDists(iSlice)))));
                end
            else
                shiftPropKer = 1i * ones(ny, nx, nSlicePerStack);
                for iSlice = 1 : nSlicePerStack
                    shiftPropKer(:, :, iSlice) = fftshift(FresnelPropKernel_X(obj.lx, ...
                        obj.ly, nx, ny, wavLen, obj.sliceDists(iSlice)));
                end
            end

            if obj.useFullAberration
                otfPhase = AberrationPhaseShift_X(obj.fullAberrations, wavLen, ...
                    obj.lx, obj.ly, nx, ny);
                otf = obj.aperture .* exp(-1i * otfPhase);
            else
                params = obj.reducedAberrations;
                params.KeV = obj.keV;
                otf = obj.aperture .* ObjTransFunc_X(params, obj.lx, obj.ly, nx, ny);
            end

            % fftshift all the transmission function in place:
            tmpTransFuncs = obj.transFuncs;
            tmpTransFuncs = fftshift(tmpTransFuncs, 1);
            tmpTransFuncs = fftshift(tmpTransFuncs, 2);

            scanNx = length(obj.scanx);
            scanNy = length(obj.scany);
            outSliceInds = DepthsToSliceIndices(obj.sliceDists, obj.nStack, obj.depths);

            totalTask = scanNy * scanNx;
            wb = waitbar(0, 'start scanning');
            rsz = 1 / obj.nBin;
            for scanIy = 1 : scanNy
                for scanIx = 1 : scanNx
                    if obj.useGPU
                        tmpWave = fftshift(gpuArray(single(GenerateProbe_X(otf, ...
                            obj.scanx(scanIx), obj.scany(scanIy), obj.lx, obj.ly, ...
                            nx, ny))));
                    else
                        tmpWave = fftshift(GenerateProbe_X(otf, ...
                            obj.scanx(scanIx), obj.scany(scanIy), obj.lx, obj.ly, ...
                            nx, ny));
                    end

                    sliceCount = 0;
                    outIdx = 1;
                    for iStack = 1 : obj.nStack
                        for iSliceInStack = 1 : nSlicePerStack
                            tmpWave = tmpWave .* tmpTransFuncs(:, :, iSliceInStack);
                            tmpWave = ifft2(shiftPropKer(:, :, iSliceInStack) .* fft2(tmpWave));

                            sliceCount = sliceCount + 1;
                            if sliceCount == outSliceInds(outIdx)
                                tmpCbed = imresize(abs((ifftshift(fft2(tmpWave))).^2), rsz, "bilinear");
                                tmpCbed = tmpCbed / mean(tmpCbed, 'all');
                                filename = ['scan_x_', num2str(scanIx - 1), '_y_',...
                                    num2str(scanIy - 1), '_z=',... 
                                    num2str(obj.depths(outIdx), '%.6f'), 'A.bin'];
                                WriteBinaryFile(fullfile(obj.cbedDir, filename), ...
                                    gather(tmpCbed), 'row', 'float');
                                outIdx = outIdx + 1;
                            end
                        end
                    end

                    doneTask = (scanIy - 1) * scanNx + scanIx;
                    doneRatio = doneTask / totalTask;
                    waitbar(doneRatio, wb, [num2str(roundn(doneRatio, -3) * 100), '%']);
                end
            end

            close(wb);

        end
    end
end