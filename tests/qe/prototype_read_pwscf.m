% prototype_read_pwscf.m
clc;
clear;
close all;
%% main:
lines = readlines('tests/qe/si.1_scf.in');
nLine = numel(lines);
iLine = 1;
while iLine < nLine
    line = strtrim(lines{iLine});
    if strcmpi(line, '&CONTROL')
        pwscf.control.calculation = '';
        iLine = iLine + 1;
        tmpStr = strtrim(lines{iLine});
        while ~strcmp(tmpStr, '/')
            tmpStrList = split(tmpStr, [",", "!", "#"]);
            for ivar = 1 : numel(tmpStrList)
                if ~isempty(tmpStrList{ivar})
                    pwscf.control = ParseControlVariable(pwscf.control, tmpStrList{ivar});
                end
            end
            iLine = iLine + 1;
            tmpStr = strtrim(lines{iLine});
        end
    elseif strcmpi(line, '&SYSTEM')
        pwscf.system.ibrav = -1;
        iLine = iLine + 1;
        tmpStr = strtrim(lines{iLine});
        while ~strcmp(tmpStr, '/')
            tmpStrList = split(tmpStr, [",", "!", "#"]);
            for ivar = 1 : numel(tmpStrList)
                if ~isempty(tmpStrList{ivar})
                    pwscf.system = ParseSystemVariable(pwscf.system, tmpStrList{ivar});
                end
            end
            iLine = iLine + 1;
            tmpStr = strtrim(lines{iLine});
        end
    elseif strcmpi(line, '&ELECTRONS')
        pwscf.electrons.conv_thr = 1e-6;
        iLine = iLine + 1;
        tmpStr = strtrim(lines{iLine});
        while ~strcmp(tmpStr, '/')
            tmpStrList = split(tmpStr, [",", "!", "#"]);
            for ivar = 1 : numel(tmpStrList)
                if ~isempty(tmpStrList{ivar})
                    pwscf.electrons = ParseElectronsVariable(pwscf.electrons, tmpStrList{ivar});
                end
            end
            iLine = iLine + 1;
            tmpStr = strtrim(lines{iLine});
        end
    elseif strcmpi(line, 'ATOMIC_SPECIES')
        typeCount = 0;
        pwscf.atomic_species.types = zeros(1, pwscf.system.ntyp);
        pwscf.atomic_species.masses = zeros(1, pwscf.system.ntyp);
        pwscf.atomic_species.pseudo_pots = strings(1, pwscf.system.ntyp);
        while typeCount < pwscf.system.ntyp
            iLine = iLine + 1;
            tmpStr = strtrim(lines{iLine});
            if ~(strcmp(tmpStr(1), '#') || strcmp(tmpStr(1), '!'))
                typeCount = typeCount + 1;
                tmpStrList = split(tmpStr);
                pwscf.atomic_species.types(typeCount) = AtomTypeStrToIdx(tmpStrList(1));
                pwscf.atomic_species.masses(typeCount) = str2double(tmpStrList(2));
                pwscf.atomic_species.pseudo_pots(typeCount) = tmpStrList(3);
            end
        end
    elseif contains(line, 'ATOMIC_POSITIONS')
        if strlength(line) == 16
            pwscf.atomic_positions.option = 'alat';
        else
            tmpStr = line(17 : end);
            pwscf.atomic_positions.option = tmpStr(isletter(tmpStr)); % crystal_sg will be read as crystalsg
        end
        atomCount = 0;
        pwscf.atomic_positions.types = zeros(1, pwscf.system.nat);
        pwscf.atomic_positions.xs = zeros(1, pwscf.system.nat);
        pwscf.atomic_positions.ys = zeros(1, pwscf.system.nat);
        pwscf.atomic_positions.zs = zeros(1, pwscf.system.nat);
        pwscf.atomic_positions.if_pos1 = ones(1, pwscf.system.nat);
        pwscf.atomic_positions.if_pos2 = ones(1, pwscf.system.nat);
        pwscf.atomic_positions.if_pos3 = ones(1, pwscf.system.nat);
        while atomCount < pwscf.system.nat
            iLine = iLine + 1;
            tmpStr = strtrim(lines{iLine});
            if ~(strcmp(tmpStr(1), '#') || strcmp(tmpStr(1), '!'))
                atomCount = atomCount + 1;
                tmpStrList = split(tmpStr);
                pwscf.atomic_positions.types(atomCount) = AtomTypeStrToIdx(tmpStrList(1));
                pwscf.atomic_positions.xs(atomCount) = str2double(tmpStrList(2));
                pwscf.atomic_positions.ys(atomCount) = str2double(tmpStrList(3));
                pwscf.atomic_positions.zs(atomCount) = str2double(tmpStrList(4));
                if numel(tmpStrList) == 7
                    pwscf.atomic_positions.if_pos1(atomCount) = str2double(tmpStrList(5));
                    pwscf.atomic_positions.if_pos2(atomCount) = str2double(tmpStrList(6));
                    pwscf.atomic_positions.if_pos3(atomCount) = str2double(tmpStrList(7));
                end
            end
        end
    else
        iLine = iLine + 1;
    end
end