function [Wave] = WaveFuncRead_0(WaveFileName, Nx, Ny, SliceIdx)
%WaveFuncRead_0.m is the zeroth version of read wave function produced by
%multislice function when optional arguments are input.
%   WaveFileName -- saved wave functions;
%   Nx, Ny -- actual sampling parameters of the wave function
%   SliceIdx -- the index of wave functions, ranging from 1 to N+1, note
%       that wave function file is a 2*(N+1) by Ny*Nx matrix, N+1 denotes
%       the exit wave;

fileID = fopen(WaveFileName, 'r');
for i = 1 : SliceIdx
    WaveReal = fgetl(fileID);
    WaveImag = fgetl(fileID);
end
fclose(fileID);
WaveReal = strsplit(WaveReal);
WaveReal = WaveReal(1 : length(WaveReal) - 1);
WaveReal = str2double(WaveReal);
WaveReal = reshape(WaveReal, Ny, Nx);

WaveImag = strsplit(WaveImag);
WaveImag = WaveImag(1 : length(WaveImag) - 1);
WaveImag = str2double(WaveImag);
WaveImag = reshape(WaveImag, Ny, Nx);

Wave = WaveReal + 1i * WaveImag;

end

