function [ExitWave] = multislice(IncidentWave, WaveLength, Lx, Ly, TransFuncs, LayerDist, StackNum, WaveFileName, SavePrec)
%MULTISLICE.M performs the multislice procedure. See E. J. Kirkland
%Advanced Computing in Electron Microscopy.
%   IncidentWave;
%   Lx, Ly -- sampling parameters;
%   TransFuncs -- transmission functions;
%   LayerDist -- distance from one layer to the next layer;
%   StackNum -- number of stackings;
%   WaveFileName -- optional: save the wave function on each slice
%       including the incident wave function and exit wave function, each
%       wave function is save as two lines with input precision in the
%       destination txt file, real component for one line and imaginary
%       component for the next line, respectively. Notice that each line
%       represents a single matrix which is rearranged by columns and then
%       saved as one line;
%   SavePrec -- half optional: if the WaveFileName is input, SavePrec must
%       be given;

switch nargin
    case 7
        Wave = IncidentWave;
        LayerNum = length(LayerDist);
        for i = 1: StackNum
            for j = 1: LayerNum
                Wave = Wave .* TransFuncs(:, :, j);
                Wave = propTF_1(Wave, Lx, Ly, WaveLength, LayerDist(j));
            end
        end
        ExitWave = Wave;
    otherwise
        Prec = strcat('%.', num2str(SavePrec), 'f  \t');
        fileID = fopen(WaveFileName, 'w');
        Wave = IncidentWave;
        fprintf(fileID, Prec, real(Wave));
        fprintf(fileID, '\r\n');
        fprintf(fileID, Prec, imag(Wave));
        fprintf(fileID, '\r\n');
        LayerNum = length(LayerDist);
        for i = 1: StackNum
            for j = 1: LayerNum
                Wave = Wave .* TransFuncs(:, :, j);
                Wave = propTF_1(Wave, Lx, Ly, WaveLength, LayerDist(j));
                fprintf(fileID, Prec, real(Wave));
                fprintf(fileID, '\r\n');
                fprintf(fileID, Prec, imag(Wave));
                fprintf(fileID, '\r\n');
            end
        end
        ExitWave = Wave;
        fclose(fileID);
end

end

