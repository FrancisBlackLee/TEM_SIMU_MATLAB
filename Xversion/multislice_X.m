function [ExitWave] = multislice_X(InciWave, WavLen, Lx, Ly, TransFuncs, SliceDist, StackNum, TransFuncDir, FileExtension)
%multislice_X.m performs the multislice procedure. See E. J. Kirkland
%Advanced Computing in Electron Microscopy for more details.
%   InciWave -- incident wave;
%   Lx, Ly -- sampling side lengths;
%   TransFuncs -- 3D array including transmission functions,
%       TransFuncs(:, :, i) denotes the ith transmission function. Note
%       that for a bulk (large volume) material possibly containing too
%       many transmission functions, creating such a 3D array might cause
%       memory overflow, thus input TransFuncs = 'files', the program will
%       load transmission functions under the given path (the optional
%       input), considering that for such bulk materials, these slices do
%       not need to be looped, these files are only loaded once.
%   SliceDist -- array whose elements are distances from the identically
%       indexed slice to the next slice;
%   StackNum -- number of stackings for the slices;
%   TransFuncDir -- directory where the transmission functions are stored.
%       Also note that these TransFunc files are named in the same style,
%       name ordering is the same as the slice ordering.
%   FileExtension -- a required input if TransFuncDir is input. '*.txt' is
%       suggested.
%   Note: X denotes an experimental version.

switch nargin
    case 7
        TempWave = fftshift(InciWave);
        SliceNum = length(SliceDist);
        [Ny, Nx] = size(InciWave);
        for SliceIdx = 1 : SliceNum
            ShiftedPropKernels(:, :, SliceIdx) = fftshift(FresnelPropKernel_X(Lx, Ly, Nx, Ny, WavLen, SliceDist(SliceIdx)));
        end
        for StackIdx = 1 : StackNum
            for SliceIdx = 1 : SliceNum
                TempWave = TempWave .* fftshift(TransFuncs(:, :, SliceIdx));
                TempWave = ifft2(ShiftedPropKernels(:, :, SliceIdx) .* fft2(TempWave));
            end
        end
        ExitWave = ifftshift(TempWave);
    otherwise
        if ~isfolder(TransFuncDir)
          errorMessage = sprintf('Error: The following folder does not exist:\n%s', TransFuncDir);
          uiwait(warndlg(errorMessage));
          return;
        end
        TransFuncFiles = dir(fullfile(TransFuncDir,FileExtension));
        FileNames = {TransFuncFiles.name}';
        SortedNames = natsortfiles(FileNames);
        TempWave = fftshift(InciWave);
        [Ny, Nx] = size(InciWave);
        for FileIdx = 1 : numel(Sortednames)
            TempTransFunc = load(fullfile(TransFuncDir, SortedNames{FileIdx}));
            TempWave = TempWave .* fftshift(TempTransFunc);
            ShiftedPropKernel = fftshift(FresnelPropKernel_X(Lx, Ly, Nx, Ny, WavLen, SliceDist(FileIdx)));
            TempWave = ifft2(ShiftedPropKernel .* fft2(TempWave));
        end
        ExitWave = ifftshift(TempWave);
end

end

