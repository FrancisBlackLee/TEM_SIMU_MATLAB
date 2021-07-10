function [atomTypeIdx] = AtomTypeStrToIdx(atomTypeStr)
%AtomTypeStrToIdx() convert the string-format atom type to the index-format
%atom type, index starting from 1.
% Input:
%   atomTypeStr -- string-format atom type;
% Output:
%   atomTypeIdx -- index-format atom type;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

upperStr = upper(atomTypeStr);
switch upperStr
    case 'H'
        atomTypeIdx = 1;
    case 'HE'
        atomTypeIdx = 2;
    case 'LI'
        atomTypeIdx = 3;
    case 'BE'
        atomTypeIdx = 4;
    case 'B'
        atomTypeIdx = 5;
    case 'C'
        atomTypeIdx = 6;
    case 'N'
        atomTypeIdx = 7;
    case 'O'
        atomTypeIdx = 8;
    case 'F'
        atomTypeIdx = 9;
    case 'NE'
        atomTypeIdx = 10;
    case 'NA'
        atomTypeIdx = 11;
    case 'MG'
        atomTypeIdx = 12;
    case 'AL'
        atomTypeIdx = 13;
    case 'SI'
        atomTypeIdx = 14;
    case 'P'
        atomTypeIdx = 15;
    case 'S'
        atomTypeIdx = 16;
    case 'CL'
        atomTypeIdx = 17;
    case 'AR'
        atomTypeIdx = 18;
    case 'K'
        atomTypeIdx = 19;
    case 'CA'
        atomTypeIdx = 20;
    case 'SC'
        atomTypeIdx = 21;
    case 'TI'
        atomTypeIdx = 22;
    case 'V'
        atomTypeIdx = 23;
    case 'CR'
        atomTypeIdx = 24;
    case 'MN'
        atomTypeIdx = 25;
    case 'FE'
        atomTypeIdx = 26;
    case 'CO'
        atomTypeIdx = 27;
    case 'NI'
        atomTypeIdx = 28;
    case 'CU'
        atomTypeIdx = 29;
    case 'ZN'
        atomTypeIdx = 30;
    case 'GA'
        atomTypeIdx = 31;
    case 'GE'
        atomTypeIdx = 32;
    case 'AS'
        atomTypeIdx = 33;
    case 'SE'
        atomTypeIdx = 34;
    case 'BR'
        atomTypeIdx = 35;
    case 'KR'
        atomTypeIdx = 36;
    case 'RB'
        atomTypeIdx = 37;
    case 'SR'
        atomTypeIdx = 38;
    case 'Y'
        atomTypeIdx = 39;
    case 'ZR'
        atomTypeIdx = 40;
    case 'NB'
        atomTypeIdx = 41;
    case 'MO'
        atomTypeIdx = 42;
    case 'TC'
        atomTypeIdx = 43;
    case 'RU'
        atomTypeIdx = 44;
    case 'RH'
        atomTypeIdx = 45;
    case 'PD'
        atomTypeIdx = 46;
    case 'AG'
        atomTypeIdx = 47;
    case 'CD'
        atomTypeIdx = 48;
    case 'IN'
        atomTypeIdx = 49;
    case 'SN'
        atomTypeIdx = 50;
    case 'SB'
        atomTypeIdx = 51;
    case 'TE'
        atomTypeIdx = 52;
    case 'I'
        atomTypeIdx = 53;
    case 'XE'
        atomTypeIdx = 54;
    case 'CS'
        atomTypeIdx = 55;
    case 'BA'
        atomTypeIdx = 56;
    case 'LA'
        atomTypeIdx = 57;
    case 'CE'
        atomTypeIdx = 58;
    case 'PR'
        atomTypeIdx = 59;
    case 'ND'
        atomTypeIdx = 60;
    case 'PM'
        atomTypeIdx = 61;
    case 'SM'
        atomTypeIdx = 62;
    case 'EU'
        atomTypeIdx = 63;
    case 'GD'
        atomTypeIdx = 64;
    case 'TB'
        atomTypeIdx = 65;
    case 'DY'
        atomTypeIdx = 66;
    case 'HO'
        atomTypeIdx = 67;
    case 'ER'
        atomTypeIdx = 68;
    case 'TM'
        atomTypeIdx = 69;
    case 'YB'
        atomTypeIdx = 70;
    case 'LU'
        atomTypeIdx = 71;
    case 'HF'
        atomTypeIdx = 72;
    case 'TA'
        atomTypeIdx = 73;
    case 'W'
        atomTypeIdx = 74;
    case 'RE'
        atomTypeIdx = 75;
    case 'OS'
        atomTypeIdx = 76;
    case 'IR'
        atomTypeIdx = 77;
    case 'PT'
        atomTypeIdx = 78;
    case 'AU'
        atomTypeIdx = 79;
    case 'HG'
        atomTypeIdx = 80;
    case 'TL'
        atomTypeIdx = 81;
    case 'PB'
        atomTypeIdx = 82;
    case 'BI'
        atomTypeIdx = 83;
    case 'PO'
        atomTypeIdx = 84;
    case 'AT'
        atomTypeIdx = 85;
    case 'RN'
        atomTypeIdx = 86;
    case 'FR'
        atomTypeIdx = 87;
    case 'RA'
        atomTypeIdx = 88;
    case 'AC'
        atomTypeIdx = 89;
    case 'TH'
        atomTypeIdx = 90;
    case 'PA'
        atomTypeIdx = 91;
    case 'U'
        atomTypeIdx = 92;
    case 'NP'
        atomTypeIdx = 93;
    case 'PU'
        atomTypeIdx = 94;
    case 'AM'
        atomTypeIdx = 95;
    case 'CM'
        atomTypeIdx = 96;
    case 'BK'
        atomTypeIdx = 97;
    case 'CF'
        atomTypeIdx = 98;
    case 'ES'
        atomTypeIdx = 99;
    case 'FM'
        atomTypeIdx = 100;
    case 'MD'
        atomTypeIdx = 101;
    case 'NO'
        atomTypeIdx = 102;
    case 'LR'
        atomTypeIdx = 103;
    otherwise
        errorMessage = sprintf('Error: Invalid input atom type string');
        uiwait(warndlg(errorMessage));
        return;
end

end

